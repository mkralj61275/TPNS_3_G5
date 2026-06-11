#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulación NS2D — Canal con paredes isotérmicas
================================================
Curso de Termofluidos — 4to año Ingeniería Mecánica

Objetivo del programa
---------------------
Resolver el flujo y la transferencia de calor en un canal 2D con paredes
a temperatura constante (θ_pared = 1) y fluido entrando frío (θ_entrada = 0).

El programa entrega:
  1. Figura con DOS paneles:
       Panel superior  — perfiles de velocidad u(y) y temperatura θ(y)
                         en distintas secciones x del canal.
       Panel inferior  — temperatura de mezcla T_m(x) a lo largo del canal.
  2. Tabla por consola con T_m, dT_m/dx y u_centro en cada sección x.

Lo que el programa NO hace
--------------------------
  • No marca ni calcula la longitud de entrada.
  • No interpreta si el flujo está o no desarrollado.
  • No normaliza los perfiles.
  Eso es trabajo del alumno.

Parámetros a ajustar
--------------------
  Modificar PHYSICAL_PARAMS y NUM_PERFILES según el caso de estudio.
  Regla práctica: si los perfiles del extremo derecho del panel superior
  todavía cambian de forma, aumentar L_x.

Uso
---
  Ejecutar directamente en Spyder (F5) o desde terminal:
  $ python main_spyder3.py
"""

from cfd_models2 import NS2DCompleto, NS2DKEpsilon
from cfd_utils import run_simulation, pre_flight_check, tabla_consola

if __name__ == "__main__":
    # 1. Configuración Inicial
    MODO_LAMINAR = "LAMINAR"
    MODO_TURBULENTO = "TURBULENTO"
# ============================================================================
# Parámetros de la simulación  ← MODIFICAR AQUÍ
# ============================================================================
# --- CONFIGURACIÓN DE EXPORTACIÓN ---
    flg_dot = True  # True: (.) para Matlab/Octave | False: (,) para Excel ES
    
    # Empaquetamos elegantemente
    csv_config = {
        'sep': ',' if flg_dot else ';',
        'dec': '.' if flg_dot else ','
    }

    
    flg_mode = MODO_LAMINAR  # MODO_TURBULENTO<-- Elegir aquí

    PHYSICAL_PARAMS = {
        'V_in': 1.0,
        'L_y': 1,
        'L_x': 40.0,
        'Re': 2000.0, 
        'Pr': 1,
        'nx': 201,
        'ny': 51,
        'dt': 5e-3,
        't_final': 65.0
    }
    x_targets = [0.1, 0.5, 1.0, 5.0, 10.0, 16.0]
    num_perfiles = 9
# ###############################################################################
# #                          No modificar desde aquí                             
# ###############################################################################
    # 2. Diccionario de Modelos
    modelos = {
        MODO_LAMINAR: NS2DCompleto,
        MODO_TURBULENTO: NS2DKEpsilon
    }

    # 3. Preparación de Parámetros (Conversión Re/Pr a visc/diff)
    if PHYSICAL_PARAMS['ny'] % 2 == 0: PHYSICAL_PARAMS['ny'] += 1
    
    D_h = 2 * PHYSICAL_PARAMS['L_y']
    PHYSICAL_PARAMS['visc'] = (PHYSICAL_PARAMS['V_in'] * D_h) / PHYSICAL_PARAMS['Re']
    PHYSICAL_PARAMS['diff'] = PHYSICAL_PARAMS['visc'] / PHYSICAL_PARAMS['Pr']

    # 4. Chequeo de Lavado
    t_lavado = 1.5 * (PHYSICAL_PARAMS['L_x'] / PHYSICAL_PARAMS['V_in'])
    if PHYSICAL_PARAMS['t_final'] < t_lavado:
        print(f"\n⚠️ AVISO: t_final insuficiente para lavado completo ({t_lavado:.1f}s sugerido).")
        input("Presione ENTER para continuar...")

     #3. Pre-flight check
    if pre_flight_check(PHYSICAL_PARAMS):
         params_to_pass = PHYSICAL_PARAMS.copy()

    # 5. Instanciación y Ejecución Única
    if flg_mode in modelos:
        # Creamos la instancia del modelo pasando el diccionario de parámetros
        model = modelos[flg_mode](**PHYSICAL_PARAMS)
        
        # Llamamos a la función de control pasándole el objeto ya creado
        resultados=run_simulation(model=model, 
                   params=PHYSICAL_PARAMS, 
                   x_targets=x_targets, 
                   num_perfiles=num_perfiles,
                   export_config=csv_config)
    else:
        print(f"ERROR: Modelo {flg_mode} no disponible.")

#=============================================================================
#
#=============================================================================










# def run_simulation() -> None:
#     Re = 2*(PHYSICAL_PARAMS['V_in'] * PHYSICAL_PARAMS['L_y']) / PHYSICAL_PARAMS['visc']
#     Pr = PHYSICAL_PARAMS['visc'] / PHYSICAL_PARAMS['diff']
# # En el main_script, antes de correr la simulación:
#     ny = PHYSICAL_PARAMS['ny']
#     if ny % 2 == 0:
#         print(f"AVISO: ny={ny} es par. Se incrementará a {ny + 1} para capturar el nodo central.")
#         PHYSICAL_PARAMS['ny'] += 1
#     print()
#     print("=" * 62)
#     print("  NS2D — Canal con paredes isotérmicas")
#     print(f"  Re = {Re:.1f}   Pr = {Pr:.2f}")
#     print(f"  V_in = {PHYSICAL_PARAMS['V_in']} m/s   "
#           f"L_x = {PHYSICAL_PARAMS['L_x']} m   "
#           f"L_y = {PHYSICAL_PARAMS['L_y']} m")
#     print("=" * 62)
#     print()

#     model = NS2DCompleto(**PHYSICAL_PARAMS)
#     print(f"Pasos temporales: {model.n_steps:,}   "
#           f"(t_final = {PHYSICAL_PARAMS['t_final']} s, "
#           f"dt = {PHYSICAL_PARAMS['dt']:.0e} s)")
#     print()

#     result = model.run()
#     result.meta['V_in_used'] = PHYSICAL_PARAMS['V_in']

#     # Verificación numérica básica
#     nan_total = np.isnan(result.T).sum() + np.isnan(result.u).sum()
#     if nan_total > 0:
#         print(f"⚠️  {nan_total} NaNs detectados — posible inestabilidad numérica.")
#         print("   Reducir dt o aumentar nit puede resolver el problema.")
#     else:
#         print("Campos numéricos sin NaNs ✅")

#     # Tabla por consola
#     tabla_consola(result, n_filas=15)

#     # Figura de dos paneles
#     graficar(result, num_perfiles=NUM_PERFILES)
#     # 2. Preguntar antes de guardar
#     respuesta_guardar = input(">>> ¿Desea almacenar los resultados en archivos CSV? (s/n): ").lower()
#     if respuesta_guardar == 's':
#         # Generar sufijo único con fecha/hora y parámetros para no pisar archivos
#         timestamp = datetime.datetime.now().strftime("%H%M")
#         suffix = f"Re{int(result.Re)}_Pr{result.Pr:.1f}_{timestamp}.csv"
#         # 1. EVOLUCIÓN AXIAL (EJE CENTRAL)
#         idx_y_mid = PHYSICAL_PARAMS['ny'] // 2
#         Tm = temperatura_mezcla(result.u, result.T, result.meta['dy'])
#         dTm_dx = np.gradient(Tm, result.x)
#         # 3. Armamos el DataFrame con la "artillería" completa
#         df_axial = pd.DataFrame({
#             'x': result.x,
#             'u_central': result.u[idx_y_mid, :],
#             'p_central': result.p[idx_y_mid, :],
#             'T_central': result.T[idx_y_mid, :],
#             'T_bulk': Tm,
#             'dT_bulk_dx': dTm_dx
#         })

#         fname_axial = f"axial_{suffix}"
#         with open(fname_axial, 'w') as f:
#             f.write(f"# Modelo: {model.nombre} | Re: {result.Re:.2f} | Pr: {result.Pr:.2f}\n")
#             df_axial.to_csv(f, sep=sep_col, decimal=sep_dec, index=False)


#     columnas_perfiles = {'y': result.y}

#     for xt in x_targets:
#         if xt > result.x[-1]: continue # Evitar errores si Lx es corta
#         idx_x = np.abs(result.x - xt).argmin()
#         rx = result.x[idx_x]
#         # Guardamos u y T para cada posición
#         columnas_perfiles[f'u_x{rx:.2f}'] = result.u[:, idx_x]
#         columnas_perfiles[f'T_x{rx:.2f}'] = result.T[:, idx_x]

#         df_perfiles = pd.DataFrame(columnas_perfiles)
#         fname_perfiles = f"perfiles_{suffix}"

#         with open(fname_perfiles, 'w') as f:
#             f.write(f"# Perfiles Transversales | Re: {result.Re:.2f} | Pr: {result.Pr:.2f}\n")
#             df_perfiles.to_csv(f, sep=sep_col, decimal=sep_dec, index=False)

#     print(f"Post-procesamiento listo: {fname_axial} y {fname_perfiles}")

# if __name__ == "__main__":
#     # --- CONFIGURACIÓN DE EXPORTACIÓN ---
#     flg_dot = True  # True: (.,) para Matlab/Octave | False: (,;) para Excel ES
#     sep_col = ',' if flg_dot else ';'
#     sep_dec = '.' if flg_dot else ','

# # ============================================================================
# # Parámetros de la simulación  ← MODIFICAR AQUÍ
# # ============================================================================

# # --- Configuración de ejecución ---
#     MODO_LAMINAR = "LAMINAR"
#     MODO_TURBULENTO = "TURBULENTO"

# # Cambiá esta bandera para elegir qué simular
#     flg_mode = MODO_LAMINAR
#     PHYSICAL_PARAMS = {
#         'V_in': 1.0,     # Velocidad de entrada        [m/s] NO CAMBIAR
#         'L_y': 1.0,      # Altura del canal              [m]
#         'L_x': 20.0,     # Canal largo
#         'Re': 100.0,     # Numero de Reynolds
#         'Pr': 0.7,       # Número de Prandlt
#         'nx': 150,       # Número de nodos en la dirección axial
#         'ny': 51,        # Número de nodos en la dirección transversal
#         'dt': 1e-3,      # Paso temporal e [s]
#         't_final': 35.0  # Tiempo final de la simulasión [s]
#     }
#     NUM_PERFILES = 9   # secciones a graficar en el panel superior
#         #PERFILES TRANSVERSALES
#     x_targets = [0.1, 0.5, 1.0, 5.0, 10.0, 16]
    
# ###############################################################################
# #                          No modificar desde aquí                             
# ###############################################################################
#             # Diccionario de modelos (el "Switch" moderno de Python)
#     modelos = {
#         MODO_LAMINAR: NS2DCompleto,
#         MODO_TURBULENTO: NS2DKEpsilon
#         }
#     # 2. Correcciones de arquitectura y conversiones
#     if PHYSICAL_PARAMS['ny'] % 2 == 0: # El numero de nodos transversales debe se impar
#         PHYSICAL_PARAMS['ny'] += 1

#     # Convertimos Re y Pr a visc y diff para que NS2DCompleto los entienda
#     D_h=2*PHYSICAL_PARAMS['L_y']
#     PHYSICAL_PARAMS['visc'] = (PHYSICAL_PARAMS['V_in'] * D_h) / PHYSICAL_PARAMS['Re']
#     PHYSICAL_PARAMS['diff'] = PHYSICAL_PARAMS['visc'] / PHYSICAL_PARAMS['Pr']

#     # 3. Verificación de Lavado (Sugerencia, no bloqueo)
#     t_residencia = PHYSICAL_PARAMS['L_x'] / PHYSICAL_PARAMS['V_in']
#     t_lavado = 1.5 * t_residencia
    
#     if PHYSICAL_PARAMS['t_final'] < t_lavado:
#         print(f"\n⚠️ AVISO TÉCNICO:")
#         print(f"El t_final ({PHYSICAL_PARAMS['t_final']}s) es menor al tiempo de lavado sugerido ({t_lavado:.1f}s).")
#         print("Es muy probable que el perfil térmico no llegue a desarrollarse completamente.")
#         input("Presione ENTER para continuar de todas formas (y ver qué pasa)...")
#     # 3. Pre-flight check
#     if pre_flight_check(PHYSICAL_PARAMS):
#         # El tiempo final se calcula automáticamente para asegurar el "lavado"
#         #t_final = 1.5 * (PHYSICAL_PARAMS['L_x'] / PHYSICAL_PARAMS['V_in'])
#         # Eliminamos Re y Pr del diccionario antes de pasar a la clase 
#         # para evitar problemas si la clase no los espera en **kwargs
#         params_to_pass = PHYSICAL_PARAMS.copy()
#         # Solo si queres limpiar, pero como tenés **kwargs no es obligatorio.
#         # 4. Instanciación y corrida


#         # Instanciamos el modelo seleccionado
#         if flg_mode in modelos:
#             print(f"MODO SELECCIONADO: {flg_mode}")
#             model = modelos[flg_mode](**PHYSICAL_PARAMS)
#         else:
#             raise ValueError(f"Modelo {flg_mode} no reconocido.")
#         # Correr simulación
#         run_simulation()
        
#         # 5. Exportación de datos
#         # ... (Tu código de exportación CSV aquí)
#     else:
#         print("Simulación cancelada por el usuario.")

