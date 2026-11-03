function generar_vectores(origen)
% Genera ambos escenarios, 2020:2075. Reanuda sin repetir MAT existentes.
if nargin < 1
    origen = '/run/media/lempo/Elements/back up JG old pc avril 2026/Documents/dew Ignacio/data';
end
raiz = fileparts(fileparts(mfilename('fullpath')));
modelos = {'2.6','8.5'};
for m = 1:numel(modelos)
    for y = 2020:2075
        carpeta = fullfile(origen,['modelo' modelos{m}],num2str(y));
        assert(isfolder(carpeta),'Falta el origen: %s',carpeta);
    end
end
for m = 1:numel(modelos)
    for y = 2020:2075
        meteo = fullfile(raiz,'data',['modelo' modelos{m}],sprintf('meteo_%d.mat',y));
        vector = fullfile(raiz,'out','vectors',['model' modelos{m}],sprintf('dew_grid_%d.mat',y));
        fprintf('\n[%s] Escenario %s, anio %d\n',datestr(now()),modelos{m},y);
        if ~exist(meteo,'file')
            vectorizacion(struct('directorio_nc',fullfile(origen,['modelo' modelos{m}],num2str(y)), ...
                'modelo',['modelo' modelos{m}],'anio',y));
        end
        if ~exist(vector,'file')
            ejecutar_beysens(modelos{m},y);
        end
        % Validar cobertura, escenario y dimensiones incluso al reanudar.
        S = load(vector,'metadatos','tiempo','puntos','resultados_anio');
        esperado = (datenum(y+1,1,1)-datenum(y,1,1))*8;
        assert(S.metadatos.anio==y && strcmp(S.metadatos.beysens.modelo,modelos{m}) ...
            && numel(S.tiempo)==esperado && size(S.resultados_anio,1)==size(S.puntos,1), ...
            'Vector incompatible o incompleto: %s',vector);
        fprintf('OK %s %d: %d registros, %d puntos.\n',modelos{m},y,esperado,size(S.puntos,1));
    end
end
fprintf('COMPLETADO: 112 vectores anuales, escenarios 2.6 y 8.5, 2020:2075.\n');
end

function ejecutar_beysens(escenario,anio_objetivo)
% El script tiene su propio espacio para no modificar el bucle del lote.
beysens_modelo = escenario;
beysens_anios = anio_objetivo;
beysens;
end
