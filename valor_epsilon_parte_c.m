% =========================================================================
% delta_t.m — TPNS-3, punto 3
% Calcula el espesor de capa límite térmica delta_t(x) desde campo2D_...csv
% Criterio: T(x, delta_t) = T_centro(x) + 0.01*(T_pared - T_centro(x))
% (distancia desde la pared donde se recorrió el 99% del salto pared-núcleo)
% Compatible con MATLAB y GNU Octave.
% =========================================================================
clear all; close all;

% ---- CONFIGURAR AQUÍ -----------------------------------------------------
archivos = { 'campo2D_NS2DCompleto_Re1500_Pr0.7_2348.csv', ...
             'campo2D_NS2DCompleto_Re1500_Pr1.0_2351.csv', ...
             'campo2D_NS2DCompleto_Re1500_Pr7.0_1633.csv' };
Prs      = [0.7, 1.0, 7.0];
T_pared  = 1.0;
% ---------------------------------------------------------------------------

isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
colores  = lines(length(archivos));

figure(1); hold on; grid on;

for f = 1:length(archivos)
    % ---- Carga ----
    if isOctave
        raw = csvread(archivos{f}, 1, 0);
    else
        raw = readmatrix(archivos{f});
    end
    x_col = raw(:,1);  y_col = raw(:,2);  T_col = raw(:,6);

    x  = unique(x_col);   y = unique(y_col);
    nx = length(x);       ny = length(y);
    T  = reshape(T_col, [nx, ny])';     % (ny x nx): filas=y, columnas=x
    % Elegir hasta que valor de L_X hacer el calculo
    cols = (x <= 25.0);
    x    = x(cols);
    T    = T(:, cols);
    nx   = length(x);
    ic = floor(ny/2) + 1;               % índice del nodo central (eje)
    Ly = y(end);

    % ---- Cálculo de delta_t en cada columna x ----
    dt = nan(1, nx);
    for j = 1:nx
        Tc     = T(ic, j);
        umbral = Tc + 0.01*(T_pared - Tc);
        col    = T(1:ic, j);            % de pared inferior (y=0) al eje

        k = find(col <= umbral, 1, 'first');
        if isempty(k)
            dt(j) = Ly/2;               % la capa llegó al eje (solapada)
        elseif k == 1
            dt(j) = 0.0;
        else
            % interpolación lineal del cruce entre los nodos k-1 y k
            frac  = (col(k-1) - umbral) / (col(k-1) - col(k));
            dt(j) = y(k-1) + frac*(y(k) - y(k-1));
        end
    end

    plot(x, dt*1000, 'LineWidth', 1.8, 'Color', colores(f,:), ...
         'DisplayName', sprintf('Pr = %.1f', Prs(f)));
end

% Referencia: semialtura del canal
plot([0 max(x)], [Ly/2 Ly/2]*1000, 'k--', 'LineWidth', 1, ...
     'DisplayName', 'L_y/2 (eje)');

xlabel('x  [m]');
ylabel('\delta_t  [mm]');
title('Espesor de capa límite térmica \delta_t(x) — criterio 99% — Re = 1500');
legend('Location', 'southeast');
hold off;