from pathlib import Path
import csv

import numpy as np
from scipy.io import loadmat
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import cartopy.crs as ccrs
import rasterio
from rasterio.merge import merge
from contextlib import ExitStack


MODELOS = ("2.6", "8.5")
LIMITES = (14.8, 16.35, 43.8, 45.35)  # oeste, este, sur, norte
RAIZ = Path(__file__).resolve().parents[1]


def main():
    archivos = [(modelo, p) for modelo in MODELOS
                for p in sorted((RAIZ / "out/vectors" / f"model{modelo}").glob("dew_grid_*.mat"))]
    if not archivos:
        raise FileNotFoundError(f"No hay datos de rocío para los modelos {MODELOS}")

    oeste, este, sur, norte = LIMITES
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

    salida = RAIZ / "out/maps"
    salida.mkdir(parents=True, exist_ok=True)
    proyeccion = ccrs.PlateCarree()
    colores_tierra = LinearSegmentedColormap.from_list(
        "tierra_suave", ["#eef2df", "#e5e7cc", "#e6d8bc", "#d5c2aa"])
    # El DEM usa cota cero para el mar; se deja visible el fondo celeste.
    terreno = np.ma.masked_where(~np.isfinite(relieve[0]) | (relieve[0] <= 0), relieve[0])
    for modelo, archivo in archivos:
        anio = archivo.stem.rsplit("_", 1)[1]
        datos = loadmat(archivo, variable_names=["resultados_anio"])["resultados_anio"]
        lat, lon, dew = datos[:, :3].T
        altura = datos[:, 4]  # Altitud en metros, extraída de los DEM de data/altitude.
        validos = (np.isfinite(dew) & np.isfinite(lat) & np.isfinite(lon)
                   & (lon >= oeste) & (lon <= este) & (lat >= sur) & (lat <= norte))
        if not validos.any():
            print(f"Sin rocío válido en la región: {archivo.name}")
            continue

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
        puntos = ax.scatter(lon[validos], lat[validos], c=dew[validos], cmap="jet_r",
                            vmin=1, vmax=27,
                            s=170, edgecolors="black", linewidths=0.4,
                            transform=proyeccion, zorder=3)
        for x, y, valor in zip(lon[validos], lat[validos], altura[validos]):
            if not np.isfinite(valor):
                continue
            dx = 4 if x < oeste + 0.06 else -4 if x > este - 0.06 else 0
            dy = -14 if y > norte - 0.05 else 8
            alineacion = "left" if dx > 0 else "right" if dx < 0 else "center"
            ax.annotate(f"{valor:.0f} m", (x, y), xycoords=proyeccion._as_mpl_transform(ax),
                        xytext=(dx, dy), textcoords="offset points", ha=alineacion,
                        fontsize=8, zorder=4)
        ax.scatter(15.2314, 44.1194, marker="*", s=180, color="red",
                   edgecolors="black", transform=proyeccion, zorder=5)
        ax.text(15.25, 44.10, "Zadar", transform=proyeccion, fontweight="bold", zorder=5)
        grilla = ax.gridlines(draw_labels=True, alpha=0.3, linestyle="--")
        grilla.top_labels = grilla.right_labels = False
        fig.colorbar(puntos, ax=ax, pad=0.02, fraction=0.025, aspect=40,
                     ticks=[1, 5, 10, 15, 20, 25, 27],
                     label="Rocío acumulado [mm/año]")
        ax.set_title(f"Rocío acumulado en {anio} — modelo {modelo}\nRelieve de fondo: altitud Copernicus")
        destino = salida / f"dew-map_{modelo}_{anio}.png"
        fig.savefig(destino, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"Mapa guardado: {destino}")


if __name__ == "__main__":
    main()
