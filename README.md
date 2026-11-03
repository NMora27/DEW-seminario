El archivo beysens.m (codes/beysens.m) estima la producción de rocío en cada punto de una malla geográfica, usando datos meteorológicos cada tres horas. Después calcula los acumulados anuales y los guarda en archivos .mat.

El flujo general es:

Datos meteorológicos → punto de rocío → estimación de rocío cada 3 horas → acumulados anuales.

### 1. Configuración y archivos

Por defecto procesa el modelo '8.5' y los años 2020 a 2075:

if ~exist('beysens_modelo', 'var'), beysens_modelo = '8.5'; end
if ~exist('beysens_anios', 'var'), beysens_anios = 2020:2075; end

Puedes cambiar esos valores antes de ejecutar el script:

beysens_modelo = '8.5';
beysens_anios = 2020:2025;
beysens;

Para cada año lee:

data/modelo8.5/meteo_2020.mat

Y genera:

out/vectors/model8.5/dew_grid_2020.mat

La raíz del proyecto se obtiene a partir de la ubicación del propio script, mediante mfilename('fullpath').

### 2. Carga y comprueba los datos

M = load(entrada);
T = numel(M.tiempo);
N = size(M.puntos,1);

Aquí:

- T es la cantidad de intervalos temporales.
- N es la cantidad de puntos geográficos.
- Las matrices meteorológicas tienen tamaño T × N: cada fila representa un instante y cada columna, un punto.

Las variables principales son:

| Variable | Significado | Unidad |
|---|---|---|
| TaZ | Temperatura del aire | °C |
| HZ | Humedad relativa | % |
| WZ | Velocidad del viento | m/s |
| CloudZ | Nubosidad | oktas, de 0 a 8 |
| prZ | Precipitación del intervalo | mm/3 h |
| altitud_m | Altitud del terreno | m |

Los assert(...) detienen la ejecución si faltan campos, las dimensiones no coinciden, las unidades son incorrectas o el año no tiene la cantidad esperada de registros separados por tres horas.

Un año normal requiere 2920 registros; uno bisiesto, 2928.

### 3. Prepara la altitud

H = altitud_m'/1000;

Convierte la altitud de metros a kilómetros, que es como la utiliza la ecuación.

La transposición ' transforma el vector de altitudes de N × 1 a 1 × N. Así, MATLAB puede aplicar la altitud de cada punto a todas las filas temporales.

El script exige que la procedencia del modelo de elevación, DEM, esté registrada. No reemplaza altitudes faltantes: esos puntos conservan NaN y no tendrán una estimación válida de rocío.

### 4. Calcula el punto de rocío

El punto de rocío es la temperatura a la que el aire alcanzaría la saturación manteniendo su contenido de vapor de agua.

Se calcula con:

  aTₐ      HR
α=────+ln( ───)
  b+Tₐ     100

      bα
T_{d}=───
      a-α

Donde a = 17.27 y b = 237.7.

En el código:

alpha = a.*M.TaZ(validos_td)./(b+M.TaZ(validos_td)) ...
        + log(M.HZ(validos_td)/100);

Td(validos_td) = b.*alpha./(a-alpha);

Solo calcula Td cuando la temperatura y humedad son finitas y la humedad cumple 0 < HR ≤ 100. Los demás valores quedan como NaN.

Los operadores .*, ./ y .^ hacen operaciones elemento por elemento, porque se trabaja con matrices completas.

### 5. Estima el rocío mediante tres términos

La expresión central es:

R = term1.*term2 + term3;

term1: efecto de la altitud y del punto de rocío

term1 = 0.37 .* (1 + 0.204323.*H - 0.0238893.*H.^2 - ...
    (18.0132 - 1.04963.*H + 0.21891.*H.^2).*1e-3 .* Td);

Es una expresión empírica que combina H y Td para ajustar la estimación.

term2: efecto térmico y de la nubosidad

term2 = ((Td+273.15)/285).^4 .* (1-M.CloudZ/8);

Convierte Td a kelvin y aplica un factor de nubosidad:

- Con cielo despejado, CloudZ = 0, el factor es 1.
- Con cielo cubierto, CloudZ = 8, el factor es 0.

Así, una mayor nubosidad reduce la contribución de term1.*term2.

term3: efecto del viento y de la diferencia de temperatura

term3 = 0.045 .* (Td-M.TaZ) ...
    .* (1+100.*(1-exp(-(M.WZ/V0).^20)));

Usa V0 = 4.4 m/s como referencia. Como normalmente Td ≤ Ta, este término suele ser negativo y reduce la estimación. La potencia 20 hace que la penalización aumente con rapidez alrededor del valor de referencia del viento.

### 6. Aplica restricciones al resultado

R(R<0) = 0;
R = R/4;
R = R .* (M.prZ<0.5);

Esto significa:

1. Las estimaciones negativas se convierten en cero.
2. Se divide el resultado entre cuatro.
3. Se asigna cero cuando hay precipitación mayor o igual a 0,5 mm en el intervalo.

El divisor 4 se conserva del código Beysens anterior, según el comentario del archivo. Este script lo registra como una normalización y etiqueta la salida en mm/intervalo de 3 h; aquí no se demuestra ni se justifica físicamente ese divisor.

Después verifica viento, nubosidad, precipitación, punto de rocío y altitud. Si algún dato necesario es inválido, el resultado del intervalo queda como NaN, en lugar de cero.

### 7. Calcula los totales anuales

rocio_anual = sum(R,1)';

Suma todas las filas temporales de cada columna y entrega un total por punto geográfico.

La precipitación anual se calcula de manera semejante:

precipitacion_anual = sum(precipitacion,1)';

Las sumas son estrictas: si un punto tiene al menos un intervalo NaN, su acumulado anual también será NaN. Así se evita presentar un año incompleto como si fuera un total válido.

Luego construye:

resultados_anio = [M.puntos(:,1:2), ...
                  rocio_anual, precipitacion_anual, altitud_m];

Sus columnas son:

| Columna | Contenido |
|---|---|
| 1 | Latitud |
| 2 | Longitud |
| 3 | Rocío anual, mm |
| 4 | Precipitación anual, mm |
| 5 | Altitud, m |

### 8. Guarda los resultados

El archivo de salida incluye R, Td, los acumulados anuales, coordenadas, altitudes, tiempos y metadatos.

save(archivo,'-struct','salida','-v7');

La opción '-struct' guarda los campos de salida como variables independientes dentro del .mat. Por ejemplo, al cargarlo podrás acceder directamente a R o rocio_anual.

En términos prácticos, R(t,n) es el rocío estimado para el intervalo t en el punto n, y rocio_anual(n) es su acumulado anual, siempre que todos sus intervalos sean válidos.


----------------------------------------------------------------


La función vectorizacion.m (codes/vectorizacion.m) extrae datos meteorológicos de archivos NetCDF, selecciona una región y los organiza en matrices que después puede utilizar beysens.m. También obtiene la altitud de cada punto desde archivos de elevación .tif.

El flujo es:

NetCDF meteorológicos + TIFF de elevación → selección espacial y temporal → conversión de unidades → archivo meteo_<año>.mat.

### 1. Recibe la configuración

Aunque lo llamemos script, técnicamente es una función:

function archivo = vectorizacion(config)

Recibe una estructura config y devuelve la ruta del archivo generado.

Por ejemplo:

config = struct( ...
    'directorio_nc', '/ruta/a/netcdf', ...
    'modelo', 'modelo8.5', ...
    'anio', 2025, ...
    'region', [42 46 13 17]);

archivo = vectorizacion(config);

Las opciones predeterminadas principales son:

| Opción | Valor predeterminado | Función |
|---|---|---|
| directorio_nc | data/netcdf | Buscar archivos meteorológicos |
| directorio_dem | data/altitude | Buscar los TIFF de elevación |
| anio | 2025 | Seleccionar el año |
| region | [42 46 13 17] | Límites geográficos |
| exigir_anio_completo | true | Rechazar años incompletos |

La región se interpreta así:

[lat_min lat_max lon_min lon_max]

Por tanto, la región predeterminada corresponde a latitudes 42° a 46° norte y longitudes 13° a 17° este.

El código copia las opciones entregadas sobre los valores predeterminados y rechaza nombres de opciones desconocidos.

### 2. Determina dónde guardar el resultado

Si no indicas directorio_salida, utiliza:

data/<modelo>/meteo_<año>.mat

Por ejemplo:

data/modelo8.5/meteo_2025.mat

Si modelo está vacío, intenta deducirlo del nombre de la carpeta de entrada. Si esa carpeta se llama como un año, por ejemplo 2025, usa el nombre de su carpeta padre.

Para producir la ruta que espera beysens.m, conviene indicar explícitamente:

'modelo', 'modelo8.5'

### 3. Busca los NetCDF y las variables

ficheros = buscar_nc(c.directorio_nc);

buscar_nc recorre la carpeta de entrada y sus subcarpetas, buscando archivos .nc.

Después procesa cinco variables:

| Nombre de salida | Variable NetCDF predeterminada | Significado |
|---|---|---|
| TaZ | tas | Temperatura del aire |
| HZ | hurs | Humedad relativa |
| WZ | sfcWind | Velocidad del viento |
| CloudZ | clt | Nubosidad |
| prZ | pr | Precipitación |

Para cada una revisa todos los archivos. Si un archivo no contiene esa variable o no tiene registros del año solicitado, lo omite.

Esto permite que los datos estén repartidos entre varios archivos.

### 4. Selecciona los puntos dentro de la región

Lee las coordenadas geográficas y construye una máscara:

mascara = lat>=c.region(1) & lat<=c.region(2) & ...
          lon>=c.region(3) & lon<=c.region(4);

mascara contiene true donde el punto está dentro de los límites.

Luego:

[ix,iy] = find(mascara);
p = [lat(mascara),lon(mascara),ix,iy];

La matriz puntos queda organizada así:

| Columna | Contenido |
|---|---|
| 1 | Latitud |
| 2 | Longitud |
| 3 | Índice en rlon |
| 4 | Índice en rlat |

rlon y rlat son las dimensiones de la malla del archivo; la selección geográfica se hace con las coordenadas lat y lon.

Además, exige que todos los archivos utilizados tengan exactamente la misma malla, para que una columna siempre corresponda al mismo lugar.

### 5. Convierte la malla en una matriz de tiempo × puntos

Esta es la parte que explica el nombre vectorización.

Los datos originales tienen tres dimensiones:

tiempo × rlon × rlat

El código identifica el orden de esas dimensiones:

orden = [find(strcmp(dims,'time')), ...
         find(strcmp(dims,'rlon')), ...
         find(strcmp(dims,'rlat'))];

Lee solo el año y el rectángulo espacial necesarios. Después reorganiza y aplana:

datos = permute(datos,orden);
datos = reshape(datos,numel(sel),[]);

Finalmente conserva únicamente los puntos que pertenecen a la máscara:

datos = datos(:,sub(:));

El resultado es una matriz T × N:

- Cada fila representa un instante.
- Cada columna representa un punto geográfico.

Por ejemplo:

TaZ(10,3)

es la temperatura del décimo registro temporal en el tercer punto de puntos.

Ese mismo punto ocupa la columna 3 en humedad, viento, nubosidad y precipitación.

### 6. Maneja faltantes y convierte unidades

ncread aplica los atributos scale_factor, add_offset y _FillValue. El código también revisa missing_value y convierte esos valores a NaN.

Después llama a:

datos = convertir(datos,nombre,unidades);

Las conversiones principales son:

| Variable | Conversión |
|---|---|
| Temperatura | Kelvin → °C, restando 273.15 |
| Humedad | Fracción → porcentaje, multiplicando por 100 |
| Nubosidad | Porcentaje → oktas, multiplicando por 8/100 |
| Viento | Comprueba que esté en m/s |
| Precipitación | Tasa → cantidad correspondiente a tres horas |

Por ejemplo, para precipitación en kg m-2 s-1:

x = x*10800;

Multiplica por los 10 800 segundos de tres horas, obteniendo el equivalente en milímetros durante ese intervalo.

Si una unidad no está contemplada, detiene la ejecución para evitar interpretar mal los datos.

### 7. Ordena y sincroniza los tiempos

Los bloques extraídos de distintos archivos se concatenan y ordenan:

datos = vertcat(bloques{:});
t = vertcat(tiempos{:});
[t,orden] = sort(t);
datos = datos(orden,:);

Luego comprueba que los registros estén separados exactamente por tres horas:

assert(all(abs(diff(t)-3/24)<1e-8), ...)

Esto detecta duplicados, huecos y frecuencias incorrectas.

También asigna cada registro a un intervalo del año:

intervalos = floor((t-datenum(c.anio,1,1))*8+1e-7);

Como hay ocho intervalos diarios:

- 0: 1 de enero, de 00:00 a 03:00.
- 1: 1 de enero, de 03:00 a 06:00.
- Y así sucesivamente.

Las variables se hacen coincidir por intervalo de tres horas. No exige que sus marcas horarias originales sean idénticas si corresponden al mismo intervalo. Conserva esas marcas en tiempos_originales.

Si existen límites temporales CF, también verifica que duren tres horas y coincidan con el intervalo asignado.

Con exigir_anio_completo = true, exige 2920 registros en un año normal o 2928 en uno bisiesto.

### 8. Obtiene la altitud de cada punto

[altitud_m,archivos_dem] = leer_altitud(c.directorio_dem,puntos);

La función leer_altitud:

1. Busca archivos .tif en la carpeta de elevación.
2. Usa gdalinfo para inspeccionarlos.
3. Exige coordenadas WGS84, EPSG:4326, y una malla sin rotación.
4. Identifica los puntos cubiertos por cada TIFF.
5. Consulta su elevación con gdallocationinfo.

Toma el valor del píxel que contiene la coordenada; no interpola.

Si un punto no tiene cobertura o el píxel contiene un valor inválido, su altitud queda en NaN. Si no hay ningún TIFF, todas las altitudes quedan en NaN, pero el archivo meteorológico igualmente puede generarse.

Esto conecta con beysens.m: los puntos sin altitud no tendrán una estimación válida de rocío.

### 9. Guarda el archivo preparado

El .mat contiene:

| Campo | Contenido |
|---|---|
| TaZ, HZ, WZ, CloudZ, prZ | Matrices meteorológicas T × N |
| puntos | Coordenadas e índices de los N puntos |
| altitud_m | Altitud de cada punto |
| tiempo | Tiempos de referencia, tomados de TaZ |
| tiempos_originales | Tiempos originales de cada variable |
| limites_temporales | Límites CF, o NaN si no existen |
| intervalos_3h | Inicio y fin de cada intervalo |
| metadatos | Unidades, archivos fuente, región y convenciones |

Se guarda mediante:

save('-mat7-binary',archivo,'-struct','salida');

Los campos de salida quedan como variables independientes dentro del archivo.

### Funciones auxiliares

Las funciones al final del archivo separan tareas específicas:

- atributo: lee un atributo NetCDF o devuelve un valor predeterminado.
- leer_tiempo: convierte el tiempo NetCDF a fechas numéricas de Octave; rechaza calendarios no contemplados, como 360_day.
- leer_malla y coordenada: ordenan las coordenadas y normalizan longitudes a [-180,180).
- convertir: normaliza las unidades meteorológicas.
- leer_altitud: consulta los TIFF de elevación.
- citar: protege las rutas al construir comandos de shell.
- borrar_temporales: elimina los archivos temporales usados para consultar GDAL.
- buscar_nc: busca NetCDF recursivamente.

La diferencia principal con el archivo anterior es que vectorizacion.m prepara y valida los datos; beysens.m utiliza esos datos para calcular el rocío.


-------------------------------------------------------------

¿que metricas pueden serme utiles para determinar que puntos de todos los que calculamos son los que tienen mayor variabilidad a lo largo de los años?

Para comparar qué puntos presentan mayor variabilidad entre años, puedes trabajar con la serie de rocío anual de cada punto, es decir, los valores de rocio_anual que genera beysens.m para los distintos años.

Te recomendaría estas métricas:

| Métrica | Qué indica | Para qué te sirve |
|---|---|---|
| Desviación estándar | Cuánto se separan los valores anuales de su promedio, en mm/año | Identificar puntos con mayores fluctuaciones absolutas |
| Coeficiente de variación (CV) | Desviación estándar dividida por la media, expresada en % | Comparar la variabilidad relativa entre puntos con distinta producción |
| Rango intercuartílico (IQR) | Diferencia entre los percentiles 75 y 25 | Medir la dispersión habitual con menor influencia de años extremos |
| Cambio absoluto medio entre años consecutivos | Promedio de cuánto cambia la producción de un año al siguiente | Detectar puntos con mayor inestabilidad interanual |
| Pendiente de la tendencia | Aumento o disminución de la producción con el tiempo | Distinguir una evolución sostenida de las fluctuaciones |

### 1. Desviación estándar: variabilidad absoluta

Para un punto con producción anual R_{y}:

\[
s=\sqrt{\frac{\sum_{y=1}^{n}(R_y-\overline{R})^2}{n-1}}
\]

Una desviación estándar alta indica que la producción cambia mucho en términos de cantidad de agua.

Por ejemplo, un punto con media de 20 mm/año y desviación de 5 mm/año presenta mayores fluctuaciones absolutas que otro con media de 5 mm/año y desviación de 2 mm/año.

Usaría esta como primera métrica para ordenar los puntos.

### 2. Coeficiente de variación: variabilidad respecto a su producción

\[
CV=100\frac{s}{\overline{R}}
\]

En el ejemplo anterior:

- Punto A: CV=100(5/20)=25%.
- Punto B: CV=100(2/5)=40%.

El punto A fluctúa más en milímetros, pero el B fluctúa más en relación con su producción habitual.

Conviene mostrar el CV junto con la media: cuando la producción media está cerca de cero, el CV puede ser enorme y resultar poco informativo. Si la media es cero, no está definido.

### 3. Cambio entre años consecutivos: inestabilidad interanual

   1   n
D=───  ∑ |R_{y}-R_{y-1}|
  n-1 y=2

Esta métrica responde directamente a: «¿Cuánto cambia, en promedio, la producción de un año al siguiente?»

Es útil porque la desviación estándar no considera el orden de los años. Estas dos series tienen la misma media y desviación estándar:

A: 5, 5, 5, 15, 15, 15
B: 5, 15, 5, 15, 5, 15

Sin embargo, B cambia mucho más entre años consecutivos.

### 4. Separar tendencia y fluctuaciones

Para tus proyecciones climáticas, esta distinción es especialmente útil. Un punto puede tener una desviación estándar alta porque su producción disminuye progresivamente durante décadas, aunque cambie poco de un año al siguiente.

Puedes ajustar:

R_{y}=α+β(y-y₀)+ε_{y}

Y obtener:

- Pendiente β: cambio sostenido, en mm/año por cada año transcurrido.
- Desviación estándar de los residuos ε_{y}: variabilidad alrededor de esa tendencia.

Así puedes distinguir puntos con una evolución marcada de aquellos con fluctuaciones interanuales fuertes.

### Qué usaría en tu estudio

Prepararía una tabla por punto con:

Latitud, longitud, altitud, años válidos, media, desviación estándar, CV, IQR, cambio interanual medio, pendiente y desviación residual.

Haría dos clasificaciones principales:

1. Mayor variabilidad absoluta: ordenar por desviación estándar.
2. Mayor fluctuación alrededor de la tendencia: ordenar por desviación residual.

El CV complementaría ambas para comparar puntos de alta y baja producción.

En tus archivos hay una precaución concreta: beysens.m deja el total anual en NaN si falta algún intervalo. Para comparar los puntos de forma justa, utiliza un periodo común y, preferentemente, los mismos años válidos. Para el cambio interanual, calcula diferencias únicamente entre años realmente consecutivos. Analiza también cada escenario —por ejemplo, 2.6 y 8.5— por separado.
