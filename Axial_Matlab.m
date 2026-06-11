% =========================================================================
% campos2D.m — Carga y graficación del archivo campo2D (Anexo NS-C corregido)
% Compatible con MATLAB y GNU Octave.
% =========================================================================
clear all; close all;

isOctave  = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
sDataFile = 'campo2D_NS2DCompleto_Re1500_Pr0.7_2348.csv';

if isOctave
    % ---------------- Rama Octave ----------------
    raw_data = csvread(sDataFile, 1, 0);
    x_col = raw_data(:, 1);
    y_col = raw_data(:, 2);
    u_col = raw_data(:, 3);
    v_col = raw_data(:, 4);
    p_col = raw_data(:, 5);
    T_col = raw_data(:, 6);
    disp('Datos cargados en Octave mediante csvread.');

    nx = length(unique(x_col));
    ny = length(unique(y_col));
    X = reshape(x_col, [nx, ny])';
    Y = reshape(y_col, [nx, ny])';
    U = reshape(u_col, [nx, ny])';
    V = reshape(v_col, [nx, ny])';
    P = reshape(p_col, [nx, ny])';
    T = reshape(T_col, [nx, ny])';
else
    % ---------------- Rama MATLAB ----------------
    try
        T_table = readtable(sDataFile);
        x_col = T_table.x; y_col = T_table.y; u_col = T_table.u;
        v_col = T_table.v; p_col = T_table.p; T_col = T_table.T;
        disp('Datos cargados en MATLAB mediante readtable.');
    catch
        raw = csvread(sDataFile, 1, 0);
        x_col = raw(:,1); y_col = raw(:,2); u_col = raw(:,3);
        v_col = raw(:,4); p_col = raw(:,5); T_col = raw(:,6);
        disp('Datos cargados en MATLAB mediante csvread.');
    end

    % Reshape también en MATLAB: contourf necesita matrices 2D
    nx = length(unique(x_col));
    ny = length(unique(y_col));
    X = reshape(x_col, [nx, ny])';
    Y = reshape(y_col, [nx, ny])';
    U = reshape(u_col, [nx, ny])';
    V = reshape(v_col, [nx, ny])';
    P = reshape(p_col, [nx, ny])';
    T = reshape(T_col, [nx, ny])';
end

% ---------------- Recorte del dominio ----------------
% Nos quedamos solo con la zona de desarrollo: x <= L_corte
L_corte = 25.0;                 % [m] longitud a graficar
cols = X(1, :) <= L_corte;      % columnas de la malla con x <= L_corte
X = X(:, cols);  Y = Y(:, cols);
U = U(:, cols);  V = V(:, cols);
P = P(:, cols);  T = T(:, cols);

% ---------------- Gráficos ----------------
figure(1);

subplot(2,1,1);
contourf(X, Y, U, 20, 'LineColor', 'none');
colormap(gca, hot);
colorbar;
title('Campo de Velocidad Axial u(x,y)');
xlabel('x [m]'); ylabel('y [m]');
axis tight;

subplot(2,1,2);
contourf(X, Y, T, 20, 'LineStyle', 'none');
colormap(gca, jet);
colorbar;
hold on
skip  = 5;
XS = X(1:skip:end, 1:skip:end);
YS = Y(1:skip:end, 1:skip:end);
US = U(1:skip:end, 1:skip:end);
VS = V(1:skip:end, 1:skip:end);
scale = 0.05;   % achicado: con Ly=0.1 m, scale=0.4 tapa el canal
US = US*scale;
VS = VS*scale;
h = quiver(XS, YS, US, VS, 0, 'k');
set(h, 'LineWidth', 1.2, 'MaxHeadSize', 0.05);
title('Campo de temperatura y vectores de Velocidad');
xlabel('x [m]'); ylabel('y [m]');
axis tight;
hold off