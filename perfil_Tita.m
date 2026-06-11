% =========================================================================
% perfiles_theta.m — TPNS-3, punto 2
% Lee campo2D_...csv (columnas: x, y, u, v, p, T) y grafica los perfiles
% transversales de temperatura adimensional Theta(y) en varias estaciones x.
% Compatible con MATLAB y GNU Octave.
% =========================================================================
clear all; close all;

% ---- CONFIGURAR AQUÍ -----------------------------------------------------
sDataFile = 'campo2D_NS2DCompleto_Re1500_Pr7.0_1633.csv';
x_est     = [0.2 1 2 5 10 15 25 30 35 40];     % estaciones x a graficar [m]
sTitulo   = 'Re = 1500,  Pr = 0.7';
% ---------------------------------------------------------------------------

isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);

% ---- Carga del CSV --------------------------------------------------------
if isOctave
    raw = csvread(sDataFile, 1, 0);
    x_col = raw(:,1);  y_col = raw(:,2);  T_col = raw(:,6);
    disp('Datos cargados en Octave mediante csvread.');
else
    try
        Ttab  = readtable(sDataFile);
        x_col = Ttab.x;  y_col = Ttab.y;  T_col = Ttab.T;
        disp('Datos cargados en MATLAB mediante readtable.');
    catch
        raw = csvread(sDataFile, 1, 0);
        x_col = raw(:,1);  y_col = raw(:,2);  T_col = raw(:,6);
        disp('Datos cargados en MATLAB mediante csvread.');
    end
end

% ---- Reconstrucción de la matriz Theta (filas = y, columnas = x) ----------
x  = unique(x_col);          % vector axial       (nx x 1)
y  = unique(y_col);          % vector transversal (ny x 1)
nx = length(x);
ny = length(y);
T  = reshape(T_col, [nx, ny])';   % Theta(fila=y, columna=x)

% ---- Gráfico de perfiles ---------------------------------------------------
figure(1); hold on; grid on;
colores = jet(length(x_est));
leyenda = cell(1, length(x_est));

for k = 1:length(x_est)
    [~, j] = min(abs(x - x_est(k)));          % columna más cercana
    plot(T(:, j), y*1000, 'LineWidth', 1.6, 'Color', colores(k,:));
    leyenda{k} = sprintf('x = %.1f m', x(j));
end

xlabel('\Theta  [-]');
ylabel('y  [mm]');
title(['Evolución del perfil \Theta(y) — ', sTitulo]);
legend(leyenda, 'Location', 'eastoutside');
xlim([-0.02 1.02]);
hold off;