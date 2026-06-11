#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from __future__ import annotations
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
#from cfd_utils import temperatura_mezcla
import numpy as np

# ===========================================================================
# Funciones Auxiliares
# ===========================================================================

# def temperatura_mezcla(u: np.ndarray, T: np.ndarray, dy: float) -> np.ndarray:
#     """Calcula la temperatura de mezcla (bulk) en cada sección x."""
#     int_uT = np.trapz(u * T, dx=dy, axis=0)
#     int_u  = np.trapz(u,     dx=dy, axis=0)
#     with np.errstate(invalid='ignore'):
#         Tm = np.where(int_u > 1e-12, int_uT / int_u, 0.0)
#     return Tm

# ===========================================================================
# Resultado estándar
# ===========================================================================

@dataclass
class CFDResult:
    u: np.ndarray
    v: np.ndarray
    p: np.ndarray
    T: np.ndarray
    y: np.ndarray
    x: np.ndarray
    Re: float
    Pr: float
    meta: dict = field(default_factory=dict)

# ===========================================================================
# Clase base
# ===========================================================================

class CFDModel(ABC):
    def __init__(self, V_in, L_y, L_x, visc, diff, t_final, dt):
        self.V_in = V_in
        self.L_y = L_y
        self.L_x = L_x
        self.visc = visc
        self.diff = diff
        self.t_final = t_final
        self.dt = dt

    nombre: str = "Modelo base"
    descripcion: str = ""

    @abstractmethod
    def run(self) -> CFDResult:
        pass

# ===========================================================================
# Modelo 1 — Navier-Stokes 2D Completo (Laminar)
# ===========================================================================

class NS2DCompleto(CFDModel):
    def __init__(self, V_in=1.0, L_y=1.0, L_x=3.0, visc=0.01, diff=0.01, 
                 t_final=2.0, dt=5e-4, nx=120, ny=51, nit=50, rho=1.0, **kwargs):
        
        if ny % 2 == 0: ny += 1 # Asegurar eje central
        super().__init__(V_in, L_y, L_x, visc, diff, t_final, dt)
        self.nx, self.ny, self.nit, self.rho = nx, ny, nit, rho
        self.n_steps = int(t_final / dt)
        self.nombre = "Navier-Stokes 2D Completo"

    def run(self) -> CFDResult:
        nx, ny = self.nx, self.ny
        dx, dy = self.L_x / (nx - 1), self.L_y / (ny - 1)

        u = np.full((ny, nx), self.V_in)
        v = np.zeros((ny, nx))
        p = np.zeros((ny, nx))
        T = np.zeros((ny, nx))
        T[0, :], T[-1, :] = 1.0, 1.0 # Paredes calientes
        n_caracteres = 20
        
        for step in range(self.n_steps):
            un, vn, Tn = u.copy(), v.copy(), T.copy()

            # Poisson
            b = np.zeros((ny, nx))
            b[1:-1, 1:-1] = (self.rho / self.dt) * (
                (un[1:-1, 2:] - un[1:-1, :-2]) / (2 * dx) + (vn[2:, 1:-1] - vn[:-2, 1:-1]) / (2 * dy)
            )
            for _ in range(self.nit):
                pn = p.copy()
                p[1:-1, 1:-1] = ((pn[1:-1, 2:] + pn[1:-1, :-2]) * dy**2 +
                                 (pn[2:, 1:-1] + pn[:-2, 1:-1]) * dx**2 -
                                 b[1:-1, 1:-1] * dx**2 * dy**2) / (2 * (dx**2 + dy**2))
                p[:, -1] = 0; p[:, 0] = p[:, 1]; p[0, :] = p[1, :]; p[-1, :] = p[-2, :]

            # Momento U (con Upwind)
            adv_u_y = np.where(vn[1:-1, 1:-1] > 0, (un[1:-1, 1:-1]-un[:-2, 1:-1])/dy, (un[2:, 1:-1]-un[1:-1, 1:-1])/dy)
            u[1:-1, 1:-1] = (un[1:-1, 1:-1] - un[1:-1, 1:-1]*(self.dt/dx)*(un[1:-1, 1:-1]-un[1:-1, :-2]) 
                             - vn[1:-1, 1:-1]*self.dt*adv_u_y - (self.dt/(2*self.rho*dx))*(p[1:-1, 2:]-p[1:-1, :-2])
                             + self.visc*self.dt/dx**2 * (un[1:-1, 2:] - 2*un[1:-1, 1:-1] + un[1:-1, :-2])
                             + self.visc*self.dt/dy**2 * (un[2:, 1:-1] - 2*un[1:-1, 1:-1] + un[:-2, 1:-1]))

            # Momento V (con Upwind)
            adv_v_y = np.where(vn[1:-1, 1:-1] > 0, (vn[1:-1, 1:-1]-vn[:-2, 1:-1])/dy, (vn[2:, 1:-1]-vn[1:-1, 1:-1])/dy)
            v[1:-1, 1:-1] = (vn[1:-1, 1:-1] - un[1:-1, 1:-1]*(self.dt/dx)*(vn[1:-1, 1:-1]-vn[1:-1, :-2]) 
                             - vn[1:-1, 1:-1]*self.dt*adv_v_y - (self.dt/(2*self.rho*dy))*(p[2:, 1:-1]-p[:-2, 1:-1])
                             + self.visc*self.dt/dx**2 * (vn[1:-1, 2:] - 2*vn[1:-1, 1:-1] + vn[1:-1, :-2])
                             + self.visc*self.dt/dy**2 * (vn[2:, 1:-1] - 2*vn[1:-1, 1:-1] + vn[:-2, 1:-1]))

            # Energía T (Upwind en X e Y) ---
            # Advección en X (suponiendo flujo mayormente hacia la derecha u > 0)
            adv_T_x = (Tn[1:-1, 1:-1] - Tn[1:-1, :-2]) / dx
            
            # Advección en Y (Upwind dinámico según el signo de v)
            adv_T_y = np.where(vn[1:-1, 1:-1] > 0, 
                               (Tn[1:-1, 1:-1] - Tn[:-2, 1:-1]) / dy, 
                               (Tn[2:, 1:-1] - Tn[1:-1, 1:-1]) / dy)
            
            # Difusión (Centrada)
            diff_T = self.diff * (
                (Tn[1:-1, 2:] - 2*Tn[1:-1, 1:-1] + Tn[1:-1, :-2]) / dx**2 +
                (Tn[2:, 1:-1] - 2*Tn[1:-1, 1:-1] + Tn[:-2, 1:-1]) / dy**2
            )

            # Update
            T[1:-1, 1:-1] = (Tn[1:-1, 1:-1] 
                             - un[1:-1, 1:-1] * self.dt * adv_T_x 
                             - vn[1:-1, 1:-1] * self.dt * adv_T_y 
                             + self.dt * diff_T)
            
            # BCs
            u[:, 0] = self.V_in; v[:, 0] = 0; T[:, 0] = 0
            u[:, -1] = u[:, -2]; v[:, -1] = 0; T[:, -1] = T[:, -2]; p[:, -1] = 0
            u[0, :], u[-1, :] = 0, 0; v[0, :], v[-1, :] = 0, 0
            T[0, :], T[-1, :] = 1.0, 1.0
            # --- Bucle Temporal ---
            # --- Indicador de Progreso (Barra de 20 caracteres) ---
            
            progreso = (step + 1) / self.n_steps
            n_estrellas = int(progreso * n_caracteres)
            n_guiones = n_caracteres - n_estrellas
            
            # Imprimimos la barra. 
            # '\r' hace que el cursor vuelva al inicio de la línea en la consola
            # 'end=""' evita que salte a una línea nueva
            barra = "#" + "*" * n_estrellas + "-" * n_guiones + "#"
            print(f"\rTrabajando: {barra} {progreso*100:5.1f}%", end="")
            
        # Al terminar el bucle, un print vacío para que el siguiente texto no se pegue
        print("\nSimulación finalizada.")
        return CFDResult(u=u, v=v, p=p, T=T, 
                         y=np.linspace(0, self.L_y, ny), 
                         x=np.linspace(0, self.L_x, nx),
                         Re=(self.V_in * 2 * self.L_y)/self.visc, 
                         Pr=self.visc/self.diff,
                         meta={'dy': dy, 'nx': nx, 'ny': ny}
                         )

# ===========================================================================
# Modelo 2 — NS 2D + Turbulencia k-ε
# ===========================================================================

class NS2DKEpsilon(CFDModel):
    C_MU = 0.09; C1E = 1.44; C2E = 1.92
    SIGMA_K = 1.0; SIGMA_E = 1.3; SIGMA_T = 0.9

    def __init__(self, V_in=1.0, L_y=1.0, 
                 L_x=3.0, visc=1e-5, 
                 diff=1e-5,
                 t_final=5.0, dt=1e-4, 
                 nx=120, ny=51, nit=50, **kwargs):
        if ny % 2 == 0: ny += 1
        super().__init__(V_in, L_y, L_x, visc, diff, t_final, dt)
        self.nx, self.ny, self.nit = nx, ny, nit
        self.n_steps = int(t_final / dt)
        self.k_init = kwargs.get('k_init', 1.5 * (0.05 * V_in)**2)
        self.eps_init = kwargs.get('eps_init', 0.09**0.75 * self.k_init**1.5 / (0.07 * L_y))
        self.k_min, self.eps_min = 1e-12, 1e-12

    def run(self) -> CFDResult:
        nx, ny = self.nx, self.ny
        dx, dy = self.L_x / (nx - 1), self.L_y / (ny - 1)
        u = np.full((ny, nx), self.V_in); 
        v = np.zeros((ny, nx)); 
        p = np.zeros((ny, nx))
        T = np.zeros((ny, nx)); 
        T[0, :], T[-1, :] = 1.0, 1.0
        k = np.full((ny, nx), self.k_init); 
        eps = np.full((ny, nx), self.eps_init)
        k[0, :], k[-1, :] = 0, 0
        n_caracteres=20
        for step in range(self.n_steps):
            un, vn, Tn, kn, en = u.copy(), v.copy(), T.copy(), k.copy(), eps.copy()
            nu_t = np.clip(self.C_MU * kn**2 / np.maximum(en, self.eps_min), 0, 1000*self.visc)
            nu_eff, al_eff = self.visc + nu_t, self.diff + nu_t / self.SIGMA_T

            # Poisson
            b = (1.0 / self.dt) * ((un[1:-1, 2:]-un[1:-1, :-2])/(2*dx) + (vn[2:, 1:-1]-vn[:-2, 1:-1])/(2*dy))
            for _ in range(self.nit):
                pn = p.copy()
                p[1:-1, 1:-1] = ((pn[1:-1, 2:]+pn[1:-1, :-2])*dy**2 + (pn[2:, 1:-1]+pn[:-2, 1:-1])*dx**2 - b*dx**2*dy**2)/(2*(dx**2+dy**2))
                p[:, -1]=0; p[:, 0]=p[:, 1]; p[0, :]=p[1, :]; p[-1, :]=p[-2, :]

            ne, ae = nu_eff[1:-1, 1:-1], al_eff[1:-1, 1:-1]
            
            # Momento U, V, T (Upwind)
            for field, field_n, diff_coeff in zip([u, v, T], [un, vn, Tn], [ne, ne, ae]):
                adv_y = np.where(vn[1:-1, 1:-1]>0, (field_n[1:-1, 1:-1]-field_n[:-2, 1:-1])/dy, (field_n[2:, 1:-1]-field_n[1:-1, 1:-1])/dy)
                field[1:-1, 1:-1] = (field_n[1:-1, 1:-1] - un[1:-1, 1:-1]*self.dt*(field_n[1:-1, 1:-1]-field_n[1:-1, :-2])/dx 
                                     - vn[1:-1, 1:-1]*self.dt*adv_y 
                                     + diff_coeff*self.dt/dx**2 * (field_n[1:-1, 2:]-2*field_n[1:-1, 1:-1]+field_n[1:-1, :-2])
                                     + diff_coeff*self.dt/dy**2 * (field_n[2:, 1:-1]-2*field_n[1:-1, 1:-1]+field_n[:-2, 1:-1]))
            # Corrección presión en U y V
            u[1:-1, 1:-1] -= (self.dt/(2*dx))*(p[1:-1, 2:]-p[1:-1, :-2])
            v[1:-1, 1:-1] -= (self.dt/(2*dy))*(p[2:, 1:-1]-p[:-2, 1:-1])

            # k y epsilon
            P_k = nu_t[1:-1, 1:-1] * ((un[2:, 1:-1]-un[:-2, 1:-1])/(2*dy))**2
            for f, fn, sc in zip([k, eps], [kn, en], [self.SIGMA_K, self.SIGMA_E]):
                deff = self.visc + nu_t[1:-1, 1:-1]/sc
                adv_y = np.where(vn[1:-1, 1:-1]>0, (fn[1:-1, 1:-1]-fn[:-2, 1:-1])/dy, (fn[2:, 1:-1]-fn[1:-1, 1:-1])/dy)
                f[1:-1, 1:-1] = (fn[1:-1, 1:-1] - un[1:-1, 1:-1]*self.dt*(fn[1:-1, 1:-1]-fn[1:-1, :-2])/dx 
                                 - vn[1:-1, 1:-1]*self.dt*adv_y 
                                 + deff*self.dt/dx**2 * (fn[1:-1, 2:]-2*fn[1:-1, 1:-1]+fn[1:-1, :-2])
                                 + deff*self.dt/dy**2 * (fn[2:, 1:-1]-2*fn[1:-1, 1:-1]+fn[:-2, 1:-1]))
            
            k[1:-1, 1:-1] += self.dt * (P_k - en[1:-1, 1:-1])
            eps[1:-1, 1:-1] += self.dt * (en[1:-1, 1:-1]/np.maximum(kn[1:-1, 1:-1], self.k_min)) * (self.C1E*P_k - self.C2E*en[1:-1, 1:-1])

            # Finalizar BCs y recortes
            k = np.maximum(k, self.k_min); eps = np.maximum(eps, self.eps_min)
            for f in [u, v, T, k, eps]:
                f[:, -1] = f[:, -2]
            u[:, 0]=self.V_in; v[:, 0]=0; T[:, 0]=0; k[:, 0]=self.k_init; eps[:, 0]=self.eps_init
            u[0, :]=0; u[-1, :]=0; v[0, :]=0; v[-1, :]=0; T[0, :]=1; T[-1, :]=1; k[0, :]=0; k[-1, :]=0
            
#            if (step + 1) % (max(1, self.n_steps // 10)) == 0:
#                Tm = temperatura_mezcla(u, T, dy)[nx//2]
#                print(f"  {(step+1)/self.n_steps*100:4.1f}% | Tm(L/2)={Tm:.4f}")
            progreso = (step + 1) / self.n_steps
            n_estrellas = int(progreso * n_caracteres)
            n_guiones = n_caracteres - n_estrellas
            
            # Imprimimos la barra. 
            # '\r' hace que el cursor vuelva al inicio de la línea en la consola
            # 'end=""' evita que salte a una línea nueva
            barra = "#" + "*" * n_estrellas + "-" * n_guiones + "#"
            print(f"\rTrabajando: {barra} {progreso*100:5.1f}%", end="")
            
        # Al terminar el bucle, un print vacío para que el siguiente texto no se pegue
        print("\nSimulación finalizada.")
        return CFDResult(u=u, v=v, p=p, T=T, 
                         y=np.linspace(0, self.L_y, ny), 
                         x=np.linspace(0, self.L_x, nx),
                         Re=(self.V_in * 2 * self.L_y)/self.visc, 
                         Pr=self.visc/self.diff,
                         meta={'dy': dy, 'ny': self.ny, 
                               'dx': dx, 'nx': self.nx,
                               'nu_t':nu_t}
                         )