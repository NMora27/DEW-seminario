function archivo = vectorizacion(config)
%VECTORIZACION Extrae meteorologia CORDEX, sin Beysens ni graficos.
% archivo = vectorizacion(struct('directorio_nc', '/ruta/netcdf'))
% Requiere GNU Octave, paquete netcdf y comandos GDAL para los DEM.
% Ver codes/README.md para opciones y convenciones de salida.
if nargin < 1, config = struct(); end
pkg load netcdf;
raiz = fileparts(fileparts(mfilename('fullpath')));
c = struct('directorio_nc', fullfile(raiz,'data','netcdf'), ...
    'directorio_dem', fullfile(raiz,'data','altitude'), ...
    'directorio_salida', '', 'directorio_datos',fullfile(raiz,'data'), ...
    'modelo','', 'anio',2025, ...
    'region',[42 46 13 17], 'exigir_anio_completo',true, ...
    'variables',struct('TaZ','tas','HZ','hurs','WZ','sfcWind', ...
    'CloudZ','clt','prZ','pr'));
campos = fieldnames(config);
for k = 1:numel(campos)
    assert(isfield(c,campos{k}), 'Opcion desconocida: %s',campos{k});
    c.(campos{k}) = config.(campos{k});
end
% Por defecto: data/<modelo>/<anio>, incluso si los NetCDF son externos.
if isempty(c.directorio_salida)
    modelo = c.modelo;
    if isempty(modelo)
        origen = regexprep(c.directorio_nc,'[\\/]+$','');
        [padre,modelo] = fileparts(origen);
        if ~isempty(regexp(modelo,'^\d{4}$','once'))
            [~,modelo] = fileparts(padre);
        end
    end
    assert(~isempty(modelo) && ~any(strcmp(modelo,{'.','..'})) && ...
        isempty(regexp(modelo,'[\\/]','once')), ...
        'Indique modelo como nombre de carpeta, por ejemplo modelo8.5.');
    c.directorio_salida = fullfile(c.directorio_datos,modelo,sprintf('%d',c.anio));
end
assert(numel(c.region)==4 && c.region(1)<=c.region(2) && ...
    c.region(3)<=c.region(4), 'Region: [lat_min lat_max lon_min lon_max].');
ficheros = buscar_nc(c.directorio_nc);
assert(~isempty(ficheros),'No se encontraron NetCDF en %s.',c.directorio_nc);
nombres = {'TaZ','HZ','WZ','CloudZ','prZ'};
salida = struct(); fuentes = struct(); referencia = []; puntos = [];
tiempos_originales = struct(); limites_temporales = struct();
for v = 1:numel(nombres)
    nombre = nombres{v}; variable = c.variables.(nombre);
    bloques = {}; tiempos = {}; limites = {}; fuentes.(nombre) = {};
    for f = 1:numel(ficheros)
        ruta = fullfile(ficheros(f).folder,ficheros(f).name);
        info = ncinfo(ruta);
        if ~any(strcmp({info.Variables.Name},variable)), continue; end
        t = leer_tiempo(ruta);
        fechas = datevec(t); sel = find(fechas(:,1)==c.anio);
        if isempty(sel), continue; end
        [lat,lon] = leer_malla(ruta);
        mascara = lat>=c.region(1) & lat<=c.region(2) & ...
            lon>=c.region(3) & lon<=c.region(4);
        [ix,iy] = find(mascara);
        assert(~isempty(ix),'La region no contiene puntos CORDEX.');
        p = [lat(mascara),lon(mascara),ix,iy];
        if isempty(referencia)
            referencia = {lat,lon}; puntos = p;
        else
            assert(isequaln(lat,referencia{1}) && isequaln(lon,referencia{2}), ...
                'Mallas distintas en %s. No se permite mezclar columnas.',ruta);
        end
        vi = ncinfo(ruta,variable);
        dims = {vi.Dimensions.Name};
        orden = [find(strcmp(dims,'time')),find(strcmp(dims,'rlon')), ...
            find(strcmp(dims,'rlat'))];
        assert(numel(dims)==3 && numel(orden)==3, ...
            '%s debe tener dimensiones time, rlon y rlat.',variable);
        % Leer solo el tramo temporal y el rectangulo espacial seleccionado.
        inicio = ones(1,3); cuenta = [vi.Dimensions.Length];
        inicio(orden) = [sel(1),min(ix),min(iy)];
        cuenta(orden) = [numel(sel),max(ix)-min(ix)+1,max(iy)-min(iy)+1];
        assert(all(diff(sel)==1),'Seleccion temporal no contigua en %s.',ruta);
        datos = double(ncread(ruta,variable,inicio,cuenta));
        datos = reshape(datos,cuenta);
        % ncread aplica scale_factor/add_offset y _FillValue. missing_value
        % tambien puede existir independientemente de _FillValue.
        missing = atributo(ruta,variable,'missing_value',[]);
        escala = atributo(ruta,variable,'scale_factor',1);
        offset = atributo(ruta,variable,'add_offset',0);
        for m = 1:numel(missing)
            datos(datos==double(missing(m))*double(escala)+double(offset)) = NaN;
        end
        datos = permute(datos,orden);
        datos = reshape(datos,numel(sel),[]);
        sub = mascara(min(ix):max(ix),min(iy):max(iy));
        datos = datos(:,sub(:));
        unidades = atributo(ruta,variable,'units','');
        datos = convertir(datos,nombre,unidades);
        bloques{end+1} = datos; tiempos{end+1} = t(sel); %#ok<AGROW>
        bvar = atributo(ruta,'time','bounds','');
        if ~isempty(bvar)
            b = double(ncread(ruta,bvar));
            if size(b,1)~=2, b=b'; end
            assert(size(b,1)==2 && size(b,2)==numel(t),'Limites CF incompatibles.');
            [~,base,factor] = leer_tiempo(ruta);
            limites{end+1} = base+b(:,sel)'*factor;
        else
            limites{end+1} = nan(numel(sel),2);
        end
        fuentes.(nombre){end+1} = ruta;
    end
    assert(~isempty(bloques),'No hay datos %s para %d.',variable,c.anio);
    datos = vertcat(bloques{:}); t = vertcat(tiempos{:});
    [t,orden] = sort(t); datos = datos(orden,:);
    assert(all(abs(diff(t)-3/24)<1e-8), ...
        '%s contiene duplicados, huecos o intervalos distintos de 3 h.',variable);
    b = vertcat(limites{:}); b = b(orden,:);
    intervalos = floor((t-datenum(c.anio,1,1))*8+1e-7);
    assert(all(diff(intervalos)==1),'Intervalos de %s no consecutivos.',variable);
    if v==1
        tiempo = t; intervalos_ref = intervalos;
    else
        assert(isequal(intervalos,intervalos_ref), ...
            'Los intervalos de tres horas de %s no coinciden.',variable);
    end
    con_limites = all(isfinite(b),2);
    assert(all(abs(b(con_limites,2)-b(con_limites,1)-3/24)<1e-8), ...
        'Limites CF distintos de tres horas en %s.',variable);
    inicios = datenum(c.anio,1,1)+intervalos(con_limites)/8;
    assert(all(abs(b(con_limites,1)-inicios)<1e-8), ...
        'Limites CF no coinciden con los intervalos en %s.',variable);
    tiempos_originales.(nombre) = t;
    limites_temporales.(nombre) = b;
    salida.(nombre) = datos;
    fprintf('%s: %d tiempos x %d puntos extraidos.\n',nombre,size(datos,1),size(datos,2));
end
if c.exigir_anio_completo
    esperado = (datenum(c.anio+1,1,1)-datenum(c.anio,1,1))*8;
    assert(numel(tiempo)==esperado && ...
        tiempo(1)<datenum(c.anio,1,1)+3/24 && ...
        tiempo(end)>=datenum(c.anio+1,1,1)-3/24, ...
        'El anio %d no esta completo: se esperaban %d registros.',c.anio,esperado);
end
[altitud_m,archivos_dem] = leer_altitud(c.directorio_dem,puntos);
salida.puntos = puntos;
salida.altitud_m = altitud_m;
salida.tiempo = tiempo;
salida.tiempos_originales = tiempos_originales;
salida.limites_temporales = limites_temporales;
salida.intervalos_3h = datenum(c.anio,1,1)+[intervalos_ref,intervalos_ref+1]/8;
salida.metadatos = struct('anio',c.anio,'region',c.region, ...
    'columnas_puntos',{{'latitud','longitud','indice_rlon','indice_rlat'}}, ...
    'indices_base',1,'intervalo_horas',3,'fuentes_nc',fuentes, ...
    'fuentes_dem',{archivos_dem},'altitud_disponible',isfinite(altitud_m), ...
    'muestreo_dem','pixel que contiene la coordenada; NaN sin cobertura', ...
    'convencion_tiempo','datenum Octave; tiempo=TaZ; asociacion por intervalo de 3 h, sin interpolacion', ...
    'unidad_tiempo','dias desde 0000-01-00 (datenum)', ...
    'calendario','gregoriano', ...
    'unidades',struct('TaZ','degC','HZ','%','WZ','m/s', ...
    'CloudZ','oktas','prZ','mm/3h','altitud_m','m'));
if ~isfolder(c.directorio_salida), mkdir(c.directorio_salida); end
archivo = fullfile(c.directorio_salida,sprintf('meteo_%d.mat',c.anio));
save('-mat7-binary',archivo,'-struct','salida');
fprintf('Guardado %s: %d tiempos x %d puntos; %d altitudes disponibles.\n', ...
    archivo,numel(tiempo),size(puntos,1),sum(isfinite(altitud_m)));
end

function valor = atributo(ruta,var,nombre,defecto)
info = ncinfo(ruta,var);
if any(strcmp({info.Attributes.Name},nombre))
    valor = ncreadatt(ruta,var,nombre);
    if ischar(valor), valor = strtrim(valor(valor~=char(0))); end
else
    valor = defecto;
end
end

function [t,base,factor] = leer_tiempo(ruta)
cal = lower(atributo(ruta,'time','calendar','standard'));
assert(any(strcmp(cal,{'standard','gregorian','proleptic_gregorian'})), ...
    'Calendario %s no soportado; requiere tratamiento explicito.',cal);
u = atributo(ruta,'time','units','');
tok = regexp(u,'^(seconds|minutes|hours|days) since (.+)$','tokens','once');
assert(~isempty(tok),'Unidades temporales CF no soportadas: %s.',u);
s = regexprep(strtrim(tok{2}),'\s*(UTC|Z|\+00:00)\s*$','');
s = strrep(s,'T',' ');
partes = regexp(strtrim(s), ...
    '^(\d{4})-(\d{2})-(\d{2})(?: (\d{2}):(\d{2})(?::(\d+(?:\.\d+)?))?)?$','tokens','once');
assert(~isempty(partes),'Origen temporal no soportado: %s.',u);
valores = zeros(1,6);
for k=1:numel(partes), valores(k)=str2double(partes{k}); end
base = datenum(valores);
factor = [1 60 3600 86400];
factor = factor(strcmp(tok{1},{'seconds','minutes','hours','days'}))/86400;
t = base + double(ncread(ruta,'time'))*factor; t = t(:);
end

function [lat,lon] = leer_malla(ruta)
lat = coordenada(ruta,'lat'); lon = coordenada(ruta,'lon');
assert(isequal(size(lat),size(lon)),'Latitud y longitud incompatibles.');
lon = mod(lon+180,360)-180;
end

function x = coordenada(ruta,nombre)
i = ncinfo(ruta,nombre); d = {i.Dimensions.Name};
assert(numel(d)==2 && all(ismember({'rlon','rlat'},d)), ...
    '%s debe ser una coordenada 2D sobre rlon y rlat.',nombre);
x = double(ncread(ruta,nombre));
x = permute(x,[find(strcmp(d,'rlon')),find(strcmp(d,'rlat'))]);
end

function x = convertir(x,nombre,u)
u = lower(strtrim(u)); u = strrep(u,'**','^');
switch nombre
    case 'TaZ'
        if any(strcmp(u,{'k','kelvin'})), x=x-273.15;
        else, assert(any(strcmp(u,{'degc','degree_celsius','degrees_celsius','celsius','c','°c'})),'Unidad de temperatura desconocida: %s.',u); end
    case {'HZ','CloudZ'}
        if any(strcmp(u,{'1','fraction'})), x=x*100;
        else, assert(any(strcmp(u,{'%','percent','percentage','oktas','okta'})),'Unidad desconocida: %s.',u); end
        if strcmp(nombre,'CloudZ') && ~any(strcmp(u,{'oktas','okta'})), x=x*8/100; end
        assert(~strcmp(nombre,'HZ') || ~any(strcmp(u,{'oktas','okta'})), 'Humedad en oktas no valida.');
    case 'WZ'
        assert(any(strcmp(u,{'m s-1','m s^-1','m/s'})),'Unidad de viento desconocida: %s.',u);
    case 'prZ'
        if any(strcmp(u,{'kg m-2 s-1','kg m^-2 s^-1','mm/s'})), x=x*10800;
        elseif any(strcmp(u,{'mm/day','mm d-1','mm day-1'})), x=x/8;
        else, assert(any(strcmp(u,{'mm','mm/3h','kg m-2','kg m^-2'})),'Unidad de precipitacion desconocida: %s.',u); end
end
end

function [z,fuentes] = leer_altitud(directorio,puntos)
z = nan(size(puntos,1),1); fuentes = {};
archivos = dir(fullfile(directorio,'*.tif'));
if isempty(archivos), return; end
[estado,~] = system('gdalinfo --version');
assert(estado==0,'Se requieren gdalinfo y gdallocationinfo (GDAL).');
entrada = [tempname() '.txt']; resultado = [tempname() '.txt'];
limpieza = onCleanup(@() borrar_temporales(entrada,resultado));
for k=1:numel(archivos)
    ruta = fullfile(archivos(k).folder,archivos(k).name);
    [estado,texto] = system(['gdalinfo -json ' citar(ruta)]);
    assert(estado==0,'No se pudo inspeccionar %s.',ruta);
    info = jsondecode(texto);
    assert(~isempty(regexp(info.coordinateSystem.wkt, ...
        '(ID\["EPSG",4326\]|AUTHORITY\["EPSG","4326"\])\]$', 'once')), ...
        'El DEM debe usar coordenadas geograficas WGS84 (EPSG:4326): %s.',ruta);
    g = info.geoTransform;
    assert(g(3)==0 && g(5)==0,'DEM con rotacion no soportado: %s.',ruta);
    columnas = floor((puntos(:,2)-g(1))/g(2));
    filas = floor((puntos(:,1)-g(4))/g(6));
    indices = find(isnan(z) & columnas>=0 & columnas<info.size(1) & ...
        filas>=0 & filas<info.size(2));
    if isempty(indices), continue; end
    fid = fopen(entrada,'w');
    assert(fid>=0,'No se pudo crear archivo temporal.');
    fprintf(fid,'%.12f %.12f\n',puntos(indices,[2 1])'); fclose(fid);
    comando = ['gdallocationinfo -valonly -b 1 -wgs84 ' citar(ruta) ...
        ' < ' citar(entrada) ' > ' citar(resultado)];
    [estado,mensaje] = system(comando);
    assert(estado==0,'Error consultando DEM: %s',mensaje);
    fid = fopen(resultado,'r'); valores = nan(numel(indices),1);
    for j=1:numel(valores)
        linea = fgetl(fid);
        assert(ischar(linea),'Respuesta GDAL incompleta: %s.',ruta);
        valores(j) = str2double(strtrim(linea));
    end
    extra = fgetl(fid); fclose(fid);
    assert(~ischar(extra),'Respuesta GDAL con filas adicionales: %s.',ruta);
        if isfield(info.bands(1),'noDataValue')
        nodata = info.bands(1).noDataValue;
        if isnumeric(nodata), valores(valores==nodata)=NaN; end
    end
    valores(~isfinite(valores))=NaN;
    sel = isfinite(valores); z(indices(sel))=valores(sel);
    if any(sel), fuentes{end+1}=ruta; end
end
end

function s = citar(s)
% Comillas simples de shell, incluidas rutas con apostrofos.
s = [char(39),strrep(s,char(39),[char(39),char(34),char(39),char(34),char(39)]),char(39)];
end

function borrar_temporales(varargin)
for k=1:nargin
    if exist(varargin{k},'file'), delete(varargin{k}); end
end
end

function archivos = buscar_nc(directorio)
archivos = dir(fullfile(directorio,'*.nc'));
entradas = dir(directorio);
for k=1:numel(entradas)
    if entradas(k).isdir && ~any(strcmp(entradas(k).name,{'.','..'}))
        archivos = [archivos; buscar_nc(fullfile(directorio,entradas(k).name))];
    end
end
end
