% =========================================================================
% solapamiento.m — TPNS-3, punto 6
% Grafica la temperatura de mezcla T_m(x) y la temperatura de eje
% T_centro(x) para los tres Pr, marcando la posicion de solapamiento
% de las capas limite termicas (donde delta_t alcanza el 90% de Ly/2).
%
% Logica:
%   1) Se reconstruyen las matrices U y T del campo2D (filas=y, cols=x).
%   2) T_m(x) = integral(u*T dy) / integral(u dy)   -> regla del trapecio
%      (cuadratura robusta, como exige la Advertencia 2 de la consigna).
%   3) T_centro(x) = fila central de T.
%   4) delta_t(x) con el criterio 99% (igual que en delta_t.m) para
%      detectar el x de solapamiento.
% Compatible con MATLAB y GNU Octave.
% =========================================================================
clear all; close all;

% ---- CONFIGURAR AQUÍ -----------------------------------------------------
archivos = { 'campo2D_NS2DCompleto_Re1500_Pr0.7_2348.csv', ...
             'campo2D_NS2DCompleto_Re1500_Pr1.0_2351.csv', ...
             'campo2D_NS2DCompleto_Re1500_Pr7.0_1633.csv' };
Prs     = [0.7, 1.0, 7.0];
T_pared = 1.0;
% ---------------------------------------------------------------------------

isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);
colores  = lines(length(archivos));

figure(1); hold on; grid on;

for f = 1:length(archivos)
    % ---- Carga y reconstruccion ----
    if isOctave
        raw = csvread(archivos{f}, 1, 0);
    else
        raw = readmatrix(archivos{f});
    end
    x_col = raw(:,1); y_col = raw(:,2);
    u_col = raw(:,3); T_col = raw(:,6);

    x  = unique(x_col);  y = unique(y_col);
    nx = length(x);      ny = length(y);
    U  = reshape(u_col, [nx, ny])';
    T  = reshape(T_col, [nx, ny])';

    ic = floor(ny/2) + 1;
    Ly = y(end);

    % ---- Temperatura de mezcla por trapecios, columna a columna ----
    %   trapz integra a lo largo de la 1ra dimension (filas = y)
    Tm = trapz(y, U.*T, 1) ./ trapz(y, U, 1);   % vector 1 x nx

    % ---- Temperatura de eje ----
    Tc = T(ic, :);

    % ---- delta_t(x) para detectar solapamiento ----
    dt = nan(1, nx);
    for j = 1:nx
        umbral = Tc(j) + 0.01*(T_pared - Tc(j));
        col    = T(1:ic, j);
        k = find(col <= umbral, 1, 'first');
        if isempty(k)
            dt(j) = Ly/2;
        elseif k == 1
            dt(j) = 0;
        else
            frac  = (col(k-1) - umbral) / (col(k-1) - col(k));
            dt(j) = y(k-1) + frac*(y(k) - y(k-1));
        end
    end
    js = find(dt >= 0.9*Ly/2, 1, 'first');   % indice de solapamiento

    % ---- Trazado ----
    plot(x, Tm, '-',  'LineWidth', 1.8, 'Color', colores(f,:), ...
         'DisplayName', sprintf('Pr = %.1f: T_m', Prs(f)));
    plot(x, Tc, '--', 'LineWidth', 1.1, 'Color', colores(f,:), ...
         'HandleVisibility', 'off');
    plot(x(js), Tm(js), 'v', 'MarkerSize', 9, ...
         'MarkerFaceColor', colores(f,:), 'MarkerEdgeColor', 'k', ...
         'HandleVisibility', 'off');

    fprintf('Pr=%.1f: solapamiento en x=%.1f m | salida: Tm=%.3f, Tc=%.3f\n', ...
            Prs(f), x(js), Tm(end), Tc(end));
end

plot([0 max(x)], [1 1], 'k:', 'HandleVisibility', 'off');
xlabel('x  [m]');
ylabel('\Theta  [-]');
title('T_m(x) (llena) y T_{centro}(x) (punteada) — \nabla = solapamiento');
legend('Location', 'southeast');
ylim([0 1.08]);
hold off;