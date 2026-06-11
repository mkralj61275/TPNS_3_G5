#!/usr/bin/env python3
"""
Verificación de Continuidad — Balance de Masa
TP NS2D — Termofluidos
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ── 1. CARGAR DATOS ────────────────────────────────────────────────────────
# Reemplazá con el nombre real de tu archivo
FNAME_CAMPO = "campo2D_NS2DCompleto_Re2000_Pr1.0_1932.csv"

df = pd.read_csv(FNAME_CAMPO, comment='#')

# Reconstruir matrices 2D
nx = df['x'].nunique()
ny = df['y'].nunique()

X = df['x'].values.reshape(ny, nx)
Y = df['y'].values.reshape(ny, nx)
U = df['u'].values.reshape(ny, nx)

x_vec = X[0, :]   # vector de posiciones x
y_vec = Y[:, 0]   # vector de posiciones y
dy = y_vec[1] - y_vec[0]

# ── 2. CAUDAL EN 3 SECCIONES ───────────────────────────────────────────────
def caudal(U, x_vec, y_vec, x_target):
    """Integra u(y) en la sección x más cercana a x_target."""
    idx = np.argmin(np.abs(x_vec - x_target))
    x_real = x_vec[idx]
    u_perfil = U[:, idx]
    Q = np.trapz(u_perfil, y_vec)
    return x_real, Q, u_perfil

Lx = x_vec[-1]
secciones = [0, Lx / 2, Lx]
resultados = []

print("\n┌─────────────────────────────────────────────────────┐")
print("│         Verificación de Conservación de Masa        │")
print("├──────────────┬────────────────┬──────────────────────┤")
print("│  x objetivo  │   x real [m]   │      Q [m²/s]        │")
print("├──────────────┼────────────────┼──────────────────────┤")

for x_obj in secciones:
    x_r, Q, _ = caudal(U, x_vec, y_vec, x_obj)
    resultados.append((x_r, Q))
    print(f"│   {x_obj:6.2f} m    │   {x_r:8.4f} m   │   {Q:14.8f}       │")

print("└──────────────┴────────────────┴──────────────────────┘")

# ── 3. ERROR RELATIVO ──────────────────────────────────────────────────────
Q_entrada = resultados[0][1]
Q_salida  = resultados[-1][1]
error_rel = abs(Q_salida - Q_entrada) / Q_entrada * 100

print(f"\n  Q entrada  = {Q_entrada:.8f} m²/s")
print(f"  Q salida   = {Q_salida:.8f} m²/s")
print(f"  Error rel. = {error_rel:.4f} %")

if error_rel < 1.0:
    print("  ✅ Conservación de masa aceptable (< 1 %)")
else:
    print("  ⚠️  Error elevado — revisar nx o condiciones de borde")

# ── 4. GRÁFICO DE PERFILES ─────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 6))

colores = ['tab:blue', 'tab:orange', 'tab:green']
etiquetas = [f'x = {r[0]:.2f} m' for r in resultados]

for i, x_obj in enumerate(secciones):
    x_r, Q, u_perfil = caudal(U, x_vec, y_vec, x_obj)
    ax.plot(u_perfil, y_vec, color=colores[i],
            label=f'x = {x_r:.1f} m  (Q={Q:.4f})', linewidth=2)

ax.axvline(1.0,  color='gray', linestyle='--', linewidth=1, label='V_in = 1.0')
ax.axvline(1.5,  color='red',  linestyle=':',  linewidth=1, label='u_Poiseuille = 1.5')
ax.set_xlabel('u [m/s]')
ax.set_ylabel('y [m]')
ax.set_title('Perfiles de velocidad — Verificación de masa')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('verificacion_masa.png', dpi=150)
plt.show()
print("\n  Figura guardada: verificacion_masa.png")