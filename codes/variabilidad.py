"""Desviación estándar interanual del rocío y mapas para ambos modelos.

Ejecutar desde cualquier carpeta: python codes/variabilidad.py
Usa todos los años guardados, exige el mismo periodo en ambos modelos y
calcula la SD muestral (ddof=1). Un año faltante invalida la SD del punto.
Los CSV incluyen todos los puntos; los mapas usan la región de mapas.py.
"""

import csv
from contextlib import ExitStack
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import cartopy.crs as ccrs
import numpy as np
from scipy.io import loadmat
import rasterio
from rasterio.merge import merge

RAIZ = Path(__file__).resolve().parents[1]
MODELOS = ("2.6", "8.5")
LIMITES = (14.8, 16.35, 43.8, 45.35)  # Igual que mapas.py.
SALIDA = RAIZ / "out" / "graphics" / "variabilidad"


def calcular_sd(modelo):
    archivos = sorted((RAIZ / "out/vectors" / f"model{modelo}").glob("dew_grid_*.mat"))
    if len(archivos) < 2:
        raise ValueError(f"Modelo {modelo}: se requieren al menos dos años.")
    series = []
    coordenadas = None
    anios = []
    for archivo in archivos:
        datos = loadmat(archivo, variable_names=["resultados_anio"])["resultados_anio"]
        if datos.ndim != 2 or datos.shape[1] < 3:
            raise ValueError(f"resultados_anio incompatible: {archivo}")
        if coordenadas is None:
            coordenadas = datos[:, :2].copy()
        elif not np.array_equal(coordenadas, datos[:, :2], equal_nan=True):
            raise ValueError(f"La malla o el orden de los puntos cambió: {archivo}")
        series.append(datos[:, 2])
        anios.append(int(archivo.stem.rsplit("_", 1)[1]))
    serie = np.stack(series)
    completos = np.isfinite(serie).all(axis=0)
    sd = np.full(serie.shape[1], np.nan)
    sd[completos] = np.std(serie[:, completos], axis=0, ddof=1)
    return coordenadas, sd, np.isfinite(serie).sum(axis=0), anios


def main():
    resultados = {modelo: calcular_sd(modelo) for modelo in MODELOS}
    if resultados["2.6"][3] != resultados["8.5"][3]:
        raise ValueError("Ambos modelos deben tener los mismos años para compararlos.")
    oeste, este, sur, norte = LIMITES
    mascaras = {}
    for modelo, (coords, sd, _, _) in resultados.items():
        lat, lon = coords.T
        mascaras[modelo] = (np.isfinite(coords).all(axis=1) & np.isfinite(sd)
                            & (lon >= oeste) & (lon <= este)
                            & (lat >= sur) & (lat <= norte))
        if not mascaras[modelo].any():
            raise ValueError(f"Modelo {modelo}: no hay SD válida en la región.")
    # La misma escala permite comparar directamente los dos escenarios.
    maximo = max(sd[mascaras[m]].max() for m, (_, sd, _, _) in resultados.items())
    norma = Normalize(vmin=0, vmax=maximo if maximo > 0 else 1)
    # mapas.py usa jet_r; jet invierte su sentido: azul mínimo, rojo máximo.
    colores = "jet"
    with ExitStack() as pila:
        fuentes = [pila.enter_context(rasterio.open(p))
                   for p in sorted((RAIZ / "data/altitude").glob("*.tif"))]
        fuentes = [s for s in fuentes if s.bounds.left < este and s.bounds.right > oeste
                   and s.bounds.bottom < norte and s.bounds.top > sur]
        if not fuentes:
            raise FileNotFoundError("No hay GeoTIFF de altitud para la región")
        if any(s.crs != rasterio.crs.CRS.from_epsg(4326) for s in fuentes):
            raise ValueError("Los GeoTIFF deben estar en EPSG:4326")
        relieve, _ = merge(fuentes, bounds=(oeste, sur, este, norte),
                           res=0.002, dtype="float32", nodata=np.nan)
    terreno = np.ma.masked_where(~np.isfinite(relieve[0]) | (relieve[0] <= 0), relieve[0])
    colores_tierra = LinearSegmentedColormap.from_list(
        "tierra_suave", ["#eef2df", "#e5e7cc", "#e6d8bc", "#d5c2aa"])
    SALIDA.mkdir(parents=True, exist_ok=True)
    proyeccion = ccrs.PlateCarree()
    for modelo, (coords, sd, n_validos, anios) in resultados.items():
        with (SALIDA / f"SD_{modelo}.csv").open("w", newline="", encoding="utf-8") as archivo:
            escritor = csv.writer(archivo)
            escritor.writerow(["latitud", "longitud", "SD_mm_anio", "anios_validos", "anios_total"])
            for (lat, lon), valor, n in zip(coords, sd, n_validos):
                escritor.writerow([lat, lon, valor, n, len(anios)])
        fig, ax = plt.subplots(figsize=(12, 9), subplot_kw={"projection": proyeccion})
        ax.set_extent(LIMITES, crs=proyeccion)
        ax.set_facecolor("#e8f5fc")
        ax.imshow(terreno, extent=LIMITES, origin="upper", transform=proyeccion,
                  cmap=colores_tierra, vmin=0, vmax=1500, zorder=0)
        with (RAIZ / "data/costa.csv").open(encoding="utf-8") as costa:
            for fila in csv.reader(costa):
                if fila and fila[0].strip().startswith("LINESTRING"):
                    texto = fila[0].strip().removeprefix("LINESTRING").strip(" ()")
                    xy = np.array([p.split()[:2] for p in texto.split(",")], dtype=float)
                    ax.plot(xy[:, 0], xy[:, 1], color="black", linewidth=0.7,
                            transform=proyeccion, zorder=1)
        validos = mascaras[modelo]
        lat, lon = coords.T
        puntos = ax.scatter(lon[validos], lat[validos], c=sd[validos], cmap=colores,
                            norm=norma, s=170, edgecolors="black", linewidths=0.4,
                            transform=proyeccion, zorder=3)
        for x, y, valor in zip(lon[validos], lat[validos], sd[validos]):
            dx = 4 if x < oeste + 0.06 else -4 if x > este - 0.06 else 0
            dy = -14 if y > norte - 0.05 else 8
            ax.annotate(f"{valor:.2f}", (x, y),
                        xycoords=proyeccion._as_mpl_transform(ax),
                        xytext=(dx, dy), textcoords="offset points", fontsize=8,
                        ha="left" if dx > 0 else "right" if dx < 0 else "center", zorder=4)
        # Mostrar también los puntos cuyo historial anual está incompleto.
        sin_sd = (np.isfinite(coords).all(axis=1) & ~np.isfinite(sd)
                  & (lon >= oeste) & (lon <= este) & (lat >= sur) & (lat <= norte))
        if sin_sd.any():
            ax.scatter(lon[sin_sd], lat[sin_sd], marker="x", color="gray",
                       transform=proyeccion, label="Sin serie completa", zorder=3)
            ax.legend(loc="lower left")
        ax.scatter(15.2314, 44.1194, marker="*", s=180, color="red",
                   edgecolors="black", transform=proyeccion, zorder=5)
        ax.text(15.25, 44.10, "Zadar", transform=proyeccion, fontweight="bold", zorder=5)
        grilla = ax.gridlines(draw_labels=True, alpha=0.3, linestyle="--")
        grilla.top_labels = grilla.right_labels = False
        fig.colorbar(puntos, ax=ax, pad=0.02, fraction=0.025, aspect=40,
                     label="Desviación estándar del rocío anual [mm/año]")
        ax.set_title(f"Variabilidad interanual del rocío — modelo {modelo}\n"
                     f"{anios[0]}–{anios[-1]} · {len(anios)} años · SD muestral (n−1)")
        destino = SALIDA / f"SD_{modelo}.png"
        fig.savefig(destino, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"{destino}: {validos.sum()} puntos con SD en el mapa; "
              f"{np.isfinite(sd).sum()}/{len(sd)} con serie completa en la malla.")


if __name__ == "__main__":
    main()
