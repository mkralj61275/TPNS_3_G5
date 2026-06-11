#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio
pio.renderers.default = "browser"
import datetime
# Importamos CFDModel solo para el Type Hinting de run_simulation
from cfd_models2 import CFDModel, CFDResult 
"""
TO DO desacoplar calculo de almacenamiento.
Cambiar modo de graficación a mathplotlib
almacenar campos en vtk
"""
#===================================================================
#                              CFD
#======================================================================
def temperatura_mezcla(u: np.ndarray, T: np.ndarray, dy: float) -> np.ndarray:
        """
        Temperatura de mezcla (bulk temperature) en cada sección x.

            T_m(x) = integral(u * T, dy) / integral(u, dy)

        Integración trapezoidal columna a columna.
        """
        int_uT = np.trapz(u * T, dx=dy, axis=0)   # shape (nx,)
        int_u  = np.trapz(u,     dx=dy, axis=0)   # shape (nx,)
        with np.errstate(invalid='ignore'):
            Tm = np.where(int_u > 1e-12, int_uT / int_u, 0.0)
        return Tm
def pre_flight_check(params: dict) -> bool:
    """
    Realiza el análisis de estabilidad antes de instanciar el modelo.
    Devuelve True si el usuario acepta los riesgos o si todo está OK.
    """
    # Cálculos previos
    dx = params['L_x'] / (params['nx'] - 1)
    dy = params['L_y'] / (params['ny'] - 1)
    dt = params['dt']
    
    # Parámetros físicos necesarios para la estabilidad
    visc = 1.0 / params['Re'] # Suponiendo V_in=1 y L_y=1
    diff = visc / params['Pr']
    
    cfl_u  = params['V_in'] * dt / dx
    cfl_v  = params['V_in'] * dt / dy
    fo_x   = visc * dt / dx**2
    fo_y   = visc * dt / dy**2
    fo_T_x = diff * dt / dx**2
    fo_T_y = diff * dt / dy**2

    # Verificación de "Salud" (Criterio de advertencia: 0.5)
    checks = [cfl_u < 0.5, cfl_v < 0.5, fo_x < 0.5, fo_y < 0.5, fo_T_x < 0.5, fo_T_y < 0.5]
    all_ok = all(checks)

    print("\n┌─ Verificación de estabilidad ──────────────────────────────")
    print(f"│  dx={dx:.4f} m   dy={dy:.4f} m   dt={dt:.2e} s")
    print(f"│  CFL_u  = {cfl_u:.4f}  {'✅' if cfl_u  < 0.5 else '⚠️ ALTO' if cfl_u < 1 else '❌ INESTABLE'}")
    print(f"│  CFL_v  = {cfl_v:.4f}  {'✅' if cfl_v  < 0.5 else '⚠️ ALTO' if cfl_v < 1 else '❌ INESTABLE'}")
    print(f"│  Fo_ν x = {fo_x:.4f}  {'✅' if fo_x   < 0.5 else '❌ INESTABLE'}")
    print(f"│  Fo_ν y = {fo_y:.4f}  {'✅' if fo_y   < 0.5 else '❌ INESTABLE'}")
    print(f"│  Fo_α x = {fo_T_x:.4f}  {'✅' if fo_T_x < 0.5 else '❌ INESTABLE'}")
    print(f"│  Fo_α y = {fo_T_y:.4f}  {'✅' if fo_T_y < 0.5 else '❌ INESTABLE'}")
    print("└────────────────────────────────────────────────────────────")

    if not all_ok:
        print("¡ADVERTENCIA! Algunos parámetros superan los límites recomendados (0.5).")
    
    respuesta = input(">>> ¿Desea comenzar la simulación? (s/n): ").lower()
    return respuesta == 's'
# ============================================================================
# Tabla de postprocesamiento por consola
# ============================================================================
def tabla_consola(result, params: dict, n_filas: int = 15) -> None:
    """
    Imprime T_centro, dT_centro/dx, u_centro y du_centro/dx.
    Permite comparar directamente las longitudes de entrada hidráulica y térmica.
    """
    PHYSICAL_PARAMS=params
    x = result.x
    ny, nx = result.u.shape
    i_c = ny // 2  # Índice del nodo central (eje)

    u_c = result.u[i_c, :]
    T_c = result.T[i_c, :]
    
    # Cálculo de gradientes axiales en el eje
    du_dx = np.gradient(u_c, x)
    dT_dx = np.gradient(T_c, x)

    idx = np.linspace(0, len(x) - 1, n_filas, dtype=int)

    print()
    print("┌─────────────────────────────────────────────────────────────────────────────┐")
    print("│             Postprocesamiento — Evolución en el Eje Central                 │")
    print("├──────────┬───────────┬─────────────┬──────────────────┬─────────────────────┤")
    print("│   x [m]  │ u_centro  │  du/dx (eje)│     T_centro     │     dT/dx (eje)     │")
    print("├──────────┼───────────┼─────────────┼──────────────────┼─────────────────────┤")
    for i in idx:
        print(f"│  {x[i]:6.3f}  │  {u_c[i]:7.4f}  │  {du_dx[i]:+9.5f}  │"
              f"     {T_c[i]:8.4f}     │      {dT_dx[i]:+9.5f}      │")
    print("└──────────┴───────────┴─────────────┴──────────────────┴─────────────────────┘")
    
    u_teorico = 1.5 * PHYSICAL_PARAMS['V_in']
    print(f"  Referencia: u_inf = {PHYSICAL_PARAMS['V_in']:.2f} | "
          f"u_parab (teórico) = {u_teorico:.3f} | T_pared = 1.000")
    print()
    # ============================================================================
    # Temperatura de mezcla  T_m(x)
    # ============================================================================

# ============================================================================
# Figura de dos paneles
# ============================================================================
def graficar(result, params:dict, num_perfiles: int = 9) -> None:
    """
    Panel superior — perfiles u(y) y θ(y) en secciones seleccionadas.
                     Un punto verde marca T_m en cada sección.
    Panel inferior — temperatura de mezcla T_m(x) continua.
                     Los puntos verdes corresponden a las mismas secciones.

    Convención de colores
    ─────────────────────
      Cian          → perfil de velocidad u(y)
      Rojo          → perfil de temperatura θ(y)
      Verde (punto) → T_m en la sección
      Verde (línea) → curva T_m(x) continua
      Gris punteado → línea base de cada sección
    """
    PHYSICAL_PARAMS=params
    dy   = result.meta['dy']
    nx_r = len(result.x)
    ly   = result.y[-1]
    v_in = PHYSICAL_PARAMS['V_in']

    # ── Escala visual ──────────────────────────────────────────────────────
    # El espacio entre perfiles consecutivos se reparte: 55 % para el perfil,
    # 45 % de margen. Se escala u y T para que su valor máximo (u_Poiseuille
    # y T=1 respectivamente) ocupe ese 55 %, permitiendo comparación visual
    # directa entre ambas variables.
    espacio = result.x[-1] / (num_perfiles - 1)
    scale_u = espacio * 0.55 / (1.5 * v_in)
    scale_T = espacio * 0.55 / 1.0

    # Distribución cuadrática de secciones: más densidad en la entrada
    t_q     = np.linspace(0.0, 1.0, num_perfiles) ** 2
    indices = np.unique((t_q * (nx_r - 1)).astype(int))

    # ── T_m en todo el canal ───────────────────────────────────────────────
    Tm = temperatura_mezcla(result.u, result.T, dy)

    # ── Figura ─────────────────────────────────────────────────────────────
    fig = make_subplots(
        rows=2, cols=1,
        row_heights=[0.65, 0.35],
        shared_xaxes=True,
        vertical_spacing=0.07,
        subplot_titles=(
            "Perfiles de velocidad (cian) y temperatura (rojo)",
            "Temperatura de mezcla  T_m(x)",
        ),
    )

    # ══════════════════════════════════════════════════════════════════════
    # PANEL SUPERIOR
    # ══════════════════════════════════════════════════════════════════════
    show = {'u': True, 'T': True, 'Tm_pt': True}

    for idx in indices:
        x0 = result.x[idx]

        # Línea base de la sección
        fig.add_trace(go.Scatter(
            x=[x0, x0], y=[0.0, ly],
            mode='lines',
            line=dict(color='rgba(180,180,180,0.25)', width=0.8, dash='dot'),
            showlegend=False,
        ), row=1, col=1)

        # Perfil u(y)
        fig.add_trace(go.Scatter(
            x=x0 + result.u[:, idx] * scale_u,
            y=result.y,
            mode='lines',
            line=dict(color='cyan', width=1.8),
            name='u(y)',
            legendgroup='u',
            showlegend=show['u'],
        ), row=1, col=1)
        show['u'] = False

        # Perfil θ(y)
        fig.add_trace(go.Scatter(
            x=x0 + result.T[:, idx] * scale_T,
            y=result.y,
            mode='lines',
            line=dict(color='#ff4d4d', width=1.8),
            name='θ(y)',
            legendgroup='T',
            showlegend=show['T'],
        ), row=1, col=1)
        show['T'] = False

        # Punto T_m en la sección (dibujado a y = Ly/2 por legibilidad)
        fig.add_trace(go.Scatter(
            x=[x0 + Tm[idx] * scale_T],
            y=[ly / 2],
            mode='markers',
            marker=dict(color='#2ecc71', size=7, symbol='circle',
                        line=dict(color='white', width=1)),
            name='T_m en sección',
            legendgroup='Tm_pt',
            showlegend=show['Tm_pt'],
        ), row=1, col=1)
        show['Tm_pt'] = False

        # Etiqueta de x debajo de cada sección
        fig.add_annotation(
            x=x0, y=-0.08,
            xref='x', yref='y',
            text=f"{x0:.1f} m",
            showarrow=False,
            font=dict(size=8, color='rgba(200,200,200,0.6)'),
        )

    # ══════════════════════════════════════════════════════════════════════
    # PANEL INFERIOR 
    # ══════════════════════════════════════════════════════════════════════
    fig.add_trace(go.Scatter(
        x=result.x,
        y=Tm,
        mode='lines',
        line=dict(color='#2ecc71', width=2),
        name='T_m(x)',
    ), row=2, col=1)

    # Puntos en las mismas secciones que los perfiles
    fig.add_trace(go.Scatter(
        x=[result.x[i] for i in indices],
        y=[Tm[i]        for i in indices],
        mode='markers',
        marker=dict(color='#2ecc71', size=8, symbol='circle',
                    line=dict(color='white', width=1.2)),
        showlegend=False,
    ), row=2, col=1)

    # Referencia T_m = 1
    fig.add_hline(
        y=1.0, row=2, col=1,
        line=dict(color='rgba(255,200,80,0.45)', width=1.2, dash='dash'),
        annotation_text='T_m = 1',
        annotation_font=dict(color='rgba(255,200,80,0.8)', size=10),
        annotation_position='right',
    )

    # ── Layout ─────────────────────────────────────────────────────────────
    subtitulo = (
        f"Re = {result.Re:.1f}  ·  Pr = {result.Pr:.2f}  ·  "
        f"L_x = {result.x[-1]:.1f} m  ·  "
        f"t_final = {PHYSICAL_PARAMS['t_final']:.1f} s"
    )
    fig.update_layout(
        title=dict(
            text=f"Canal isotérmico — desarrollo hidráulico y térmico<br>"
                 f"<sup>{subtitulo}</sup>",
            x=0.5,
            font=dict(size=14),
        ),
        plot_bgcolor='#0d0d1a',
        paper_bgcolor='#0d0d1a',
        font=dict(color='rgba(220,220,220,0.9)', size=11),
        legend=dict(
            orientation='h',
            yanchor='bottom', y=1.04,
            xanchor='right',  x=1,
            bgcolor='rgba(0,0,0,0)',
        ),
    )

    fig.update_yaxes(title_text='y [m]',
                     showgrid=True, gridcolor='rgba(255,255,255,0.06)',
                     row=1, col=1)
    fig.update_yaxes(title_text='T_m  [—]',
                     range=[-0.05, 1.12],
                     showgrid=True, gridcolor='rgba(255,255,255,0.06)',
                     row=2, col=1)
    fig.update_xaxes(showgrid=True, gridcolor='rgba(255,255,255,0.06)',
                     row=1, col=1)
    fig.update_xaxes(title_text='x  [m]',
                     showgrid=True, gridcolor='rgba(255,255,255,0.06)',
                     row=2, col=1)

    fig.show()

# ============================================================================
# Punto de entrada
# ============================================================================
def run_simulation(model: CFDModel,
                   params: dict, 
                   x_targets: list, 
                   num_perfiles: int,
                   export_config: dict = None) -> None:
    """
    Recibe un modelo ya instanciado, lo corre y realiza el post-procesamiento.
    """
    PHYSICAL_PARAMS=params
    print("\n" + "=" * 62)
    print(f"  EJECUTANDO: {model.nombre}")
    print(f"  Configuración: t_final={params['t_final']}s, dt={params['dt']:.0e}s")
    print("=" * 62 + "\n")

    # 1. Ejecución del modelo
    result = model.run()
    
    # 2. Verificación de estabilidad
    nan_total = np.isnan(result.T).sum() + np.isnan(result.u).sum()
    if nan_total > 0:
        print(f"⚠️  {nan_total} NaNs detectados — posible inestabilidad numérica.")
    else:
        print("Campos numéricos sin NaNs ✅")

    # 3. Visualización básica
    tabla_consola(result, PHYSICAL_PARAMS, n_filas=15)
    graficar(result,  PHYSICAL_PARAMS, num_perfiles=num_perfiles)

    # 4. Exportación a CSV
    respuesta_guardar = input("\n>>> ¿Desea almacenar los resultados en archivos CSV? (s/n): ").lower()
    if respuesta_guardar == 's':
        if export_config is None:
            export_config = {'sep': ',', 'dec': '.'}
    
        sep_col = export_config.get('sep', ',')
        sep_dec = export_config.get('dec', '.')
        
        timestamp = datetime.datetime.now().strftime("%H%M")
        pr_val = result.Pr
        suffix = f"{model.__class__.__name__}_Re{int(result.Re)}_Pr{pr_val:.1f}_{timestamp}.csv"
        # Extraemos perfiles del eje para mayor claridad
        idx_y_mid = params['ny'] // 2
        u_c = result.u[idx_y_mid, :]
        grad_u = np.gradient(result.u, result.y, result.x)
        grad_T = np.gradient(result.T, result.y, result.x)
        T_c = result.T[idx_y_mid, :]
        p_c = result.p[idx_y_mid, :]
        
        # Temperatura de mezcla (Bulk)
        dy_eff = result.meta.get('dy', params['L_y']/(params['ny']-1))
        Tm = temperatura_mezcla(result.u, result.T, dy_eff)
        
        # Cálculo de gradientes axiales (Diferencias finitas centrales vía NumPy)
        du_dx = np.gradient(u_c, result.x)
        # El índice 0 es la pared, el índice 1 es el primer nodo de fluido
        u_pared = result.u[0, :]      # Esto es un vector de ceros [0, 0, ..., 0]
        u_fluido = result.u[1, :]     # Velocidad en la primera celda sobre la pared
        dy = result.meta['dy']        # El paso de malla vertical
        du_dy_wall = (u_fluido - u_pared) / dy
        dT_dx = np.gradient(T_c, result.x)
        T_pared = result.T[0, :]      # Esto es un vector de ceros [0, 0, ..., 0]
        T_fluido = result.T[1, :]     # Velocidad en la primera celda sobre la pared
        dT_dy_wall = (T_fluido - T_pared) / dy
        dTm_dx = np.gradient(Tm, result.x)

        # Armamos el DataFrame con todo lo necesario para el análisis posterior
        df_axial = pd.DataFrame({
            'x': result.x,
            'u_central': u_c,
            'du_central_dx': du_dx,  # Gradiente x de velocidad en el eje
            'du_dy_wall': du_dy_wall,# Gradiente y de velocidad en el eje
            'T_central': T_c,
            'dT_central_dx': dT_dx,  # Gradiente de temperatura en el eje
            'dT_dy_wall': dT_dy_wall,
            'p_central': p_c,
            'T_bulk': Tm,
            'dT_bulk_dx': dTm_dx    # Útil para el Nusselt local
        })

        fname_axial = f"axial_{suffix}"
        # Guardado con encabezado informativo estilo "legacy"
        with open(fname_axial, 'w') as f:
            f.write(f"# Modelo: {model.nombre} | Re: {result.Re:.2f} | Pr: {result.Pr:.2f}\n")
            df_axial.to_csv(f, index=False, sep=sep_col, decimal=sep_dec)
        # --- Exportación de Campo Completo (Para Quiver/Contour) ---
        # Creamos una grilla de coordenadas para cada nodo
        X_grid, Y_grid = np.meshgrid(result.x, result.y)
        
        df_campo = pd.DataFrame({
            'x': X_grid.flatten(),
            'y': Y_grid.flatten(),
            'u': result.u.flatten(),
            'v': result.v.flatten(),
            'p': result.p.flatten(),
            'T': result.T.flatten()
        })
        
        # Si es turbulento, agregamos k y epsilon
        if hasattr(result, 'k'):
            df_campo['k'] = result.k.flatten()
            df_campo['epsilon'] = result.epsilon.flatten()

        fname_campo = f"campo2D_{suffix}"
        df_campo.to_csv(fname_campo, index=False)
        
        # --- Exportación Perfiles (Tu código original con el agregado de v) ---
        columnas_perfiles = {'y': result.y}
        for xt in x_targets:
            if xt <= result.x[-1]:
                idx_x = np.abs(result.x - xt).argmin()
                rx = result.x[idx_x]
                columnas_perfiles[f'u_x{rx:.2f}'] = result.u[:, idx_x]
                columnas_perfiles[f'v_x{rx:.2f}'] = result.v[:, idx_x] # <--- Agregamos v aquí
                columnas_perfiles[f'T_x{rx:.2f}'] = result.T[:, idx_x]
        
        df_perfiles = pd.DataFrame(columnas_perfiles)
        fname_perfiles = f"perfiles_{suffix}"
        df_perfiles.to_csv(fname_perfiles, index=False)

        # # --- Exportación Perfiles ---
        # columnas_perfiles = {'y': result.y}
        # for xt in x_targets:
        #     if xt <= result.x[-1]:
        #         idx_x = np.abs(result.x - xt).argmin()
        #         rx = result.x[idx_x]
        #         columnas_perfiles[f'u_x{rx:.2f}'] = result.u[:, idx_x]
        #         columnas_perfiles[f'T_x{rx:.2f}'] = result.T[:, idx_x]
        
        # df_perfiles = pd.DataFrame(columnas_perfiles)
        # fname_perfiles = f"perfiles_{suffix}"
        # df_perfiles.to_csv(fname_perfiles, index=False)
        
        print(f"Archivos guardados: {fname_axial} y {fname_perfiles}")
        return result
