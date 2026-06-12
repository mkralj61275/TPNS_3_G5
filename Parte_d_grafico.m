% =========================================================================
% ley_potencias.m — TPNS-3, punto 5
% Calcula delta_t(x) desde un campo2D_...csv, ajusta la ley de potencias
%   delta_t = A * x^alpha   (regresion lineal en logaritmos, Anexo NS-E)
% y compara el exponente con las referencias de literatura:
%   alpha = 1/2  (capa en corriente uniforme, tipo placa plana)
%   alpha = 1/3  (capa en perfil cortante desarrollado, Leveque)
% Compatible con MATLAB y GNU Octave.
% =========================================================================
clear all; close all;

% ---- CONFIGURAR AQUÍ -----------------------------------------------------
% Usar la corrida ZOOM de la entrada (L_x=4, nx=321, ny=81, dt=2e-3, Pr=7):
sDataFile = 'campo2D_NS2DCompleto_Re1500_Pr7.0_1633.csv';
T_pared   = 1.0;
x_min     = 1.0;    % [m] inicio del ajuste (excluye la zona de entrada
                    %     contaminada por la singularidad de la esquina)
frac_max  = 0.8;    % ajustar solo donde delta_t < frac_max*(Ly/2)
                    %     (zona de crecimiento libre, sin confinamiento)
% ---------------------------------------------------------------------------

isOctave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);

% ---- Carga ----
if isOctave
    raw = csvread(sDataFile, 1, 0);
else
    raw = readmatrix(sDataFile);
end
x_col = raw(:,1);  y_col = raw(:,2);  T_col = raw(:,6);
x  = unique(x_col);   y = unique(y_col);
nx = length(x);       ny = length(y);
T  = reshape(T_col, [nx, ny])';

ic = floor(ny/2) + 1;
Ly = y(end);

% ---- delta_t(x): criterio 99% con interpolacion ----
dt = nan(1, nx);
for j = 1:nx
    Tc     = T(ic, j);
    umbral = Tc + 0.01*(T_pared - Tc);
    col    = T(1:ic, j);
    k = find(col <= umbral, 1, 'first');
    if isempty(k)
        dt(j) = Ly/2;
    elseif k == 1
        dt(j) = 0.0;
    else
        frac  = (col(k-1) - umbral) / (col(k-1) - col(k));
        dt(j) = y(k-1) + frac*(y(k) - y(k-1));
    end
end

% ---- Seleccion de la zona de ajuste ----
mask = (x(:)' > x_min) & (dt < frac_max*Ly/2) & (dt > 0);
xm   = x(mask);
dm   = dt(mask);
fprintf('Puntos de ajuste: %d  (x de %.2f a %.2f m)\n', ...
        sum(mask), xm(1), xm(end));

% ---- Ajuste por minimos cuadrados en logaritmos (Anexo NS-E) ----
log_x = log(xm(:));
log_d = log(dm(:));
Xd    = [ones(length(log_d),1), log_x];
beta  = Xd \ log_d;
A     = exp(beta(1));
alpha = beta(2);

% Coeficiente de determinacion R^2
pred = Xd * beta;
R2   = 1 - sum((log_d - pred).^2) / sum((log_d - mean(log_d)).^2);

fprintf('=== Ajuste de ley de potencias ===\n');
fprintf('delta_t = %.4f * x^%.3f      R^2 = %.4f\n', A, alpha, R2);
fprintf('Referencias: 0.500 (placa plana) | 0.333 (Leveque)\n');

% ---- Figura: lineal + log-log ----
figure(1);

subplot(1,2,1); hold on; grid on;
plot(x, dt*1000, 'r-', 'LineWidth', 1.8);
plot(xm, A*xm.^alpha*1000, 'k--', 'LineWidth', 1.5);
plot([min(x) max(x)], [Ly/2 Ly/2]*1000, 'Color', [0.5 0.5 0.5], ...
     'LineStyle', ':');
xlabel('x [m]'); ylabel('\delta_t [mm]');
title('\delta_t(x) y ajuste');
legend('Simulacion', sprintf('Ajuste: %.3f x^{%.2f}', A, alpha), ...
       'L_y/2', 'Location', 'southeast');
hold off;

subplot(1,2,2);
loglog(xm, dm*1000, 'r.', 'MarkerSize', 8); hold on; grid on;
loglog(xm, A*xm.^alpha*1000, 'k--', 'LineWidth', 1.5);
loglog(xm, dm(1)*1000*(xm/xm(1)).^0.5,   'b:', 'LineWidth', 1.5);
loglog(xm, dm(1)*1000*(xm/xm(1)).^(1/3), 'g:', 'LineWidth', 1.5);
xlabel('x [m]'); ylabel('\delta_t [mm]');
title('Escala log-log');
legend('Simulacion (zona libre)', ...
       sprintf('Ajuste: x^{%.2f}', alpha), ...
       'x^{0.5} (placa plana)', 'x^{1/3} (Leveque)', ...
       'Location', 'northwest');
hold off;