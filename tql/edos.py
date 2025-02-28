##########################################
## PARAMETROS VARIABLES
##########################################

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from tql.parametros_variables import *
from scipy.integrate import quad, solve_ivp

def pasteurizador_ode(t, y, params: Pastparams):
    """
    EDO para la temperatura del pasteurizador, Tpast(t).
    """
    Tpast = y[0]
    m_suero = params.m_suero()/3600
    f1 = params.f1
    T_in = params.get_T_in()
    T_vap = params.get_T_vap()
    cp_suero = params.get_cp_suero(T_in, f1)
    U = params.U
    A = params.get_area_past()
    m_hold = params.get_m_hold(T_in, f1)

    Q_suero = m_suero * cp_suero * (T_in - Tpast)
    Q_vapor = U * A * (T_vap - Tpast)

    dTpast_dt = (Q_suero + Q_vapor) / (m_hold * cp_suero)
    return dTpast_dt

def simulate_pasteurization(t_span, T0, params: Pastparams):
    """
    Integra la ODE del pasteurizador en [t_span].
    """
    sol = solve_ivp(
        fun=lambda t, y: pasteurizador_ode(t, y, params),
        t_span=t_span,
        y0=[T0],
        dense_output=True
    )
    return sol  # Ahora retorna la solución

def evap_ode(t, y, pastparams: Pastparams, effectparams_e1: EvapEffectParams):
    """EDO para la temperatura de la pared en el primer efecto de evaporación."""
    
    Tvap = pastparams.get_T_vap()  # Llamada correcta
    Tprod = pastparams.get_T_in()  # Este valor debería venir de una simulación previa
    Tw = y[0]
    m_prod_in = pastparams.m_suero() / 3600  # Flujo de suero en [kg/s]

    # Métodos de la instancia effectparams_e1
    h_cond = effectparams_e1.h_cond(Tvap, Tw)
    h_evap = effectparams_e1.h_evap(Tvap, Tw, m_prod_in)
    A_in = effectparams_e1.get_area_effect("in")
    A_out = effectparams_e1.get_area_effect("out")
    m_inox = effectparams_e1.m_inox()
    cp_inox = effectparams_e1.cp_inox

    if h_cond == 0:
    # Si no hay condensación, en la práctica, el calor liberado por el vapor es cero.
    # Esto permite a Tw eventualmente descender si por el lado producto (Tw -> Tprod) la pared se enfría.
        Q_cond = 0

    # Prevenir h_cond y h_evap inválidos
    if np.isnan(h_cond) or h_cond < 0:
        print(f"Advertencia: h_cond inválido ({h_cond}). Ajustando a 2200.")
        h_cond = 2200
    if np.isnan(h_evap) or h_evap < 0:
        print(f"Advertencia: h_evap inválido ({h_evap}). Ajustando a 2200.")
        h_evap = 2200

    Q_cond = h_cond * A_out * (Tvap - Tw)
    Q_evap = h_evap * A_in * (Tw - Tprod)

    # Evitar división por cero
    if m_inox * cp_inox == 0:
        print("Error: m_inox * cp_inox es 0, ajuste los parámetros.")
        return [0]

    dTdt = (Q_cond - Q_evap) / (m_inox * cp_inox)
    return [dTdt]

def simular_temp_pared(T0, t_final, pastparams, effectparams_e1):
    """Simula la ecuación diferencial para el primer efecto de evaporación."""
    sol = solve_ivp(lambda t, y: evap_ode(t, y, pastparams, effectparams_e1),
                    [0, t_final], [T0], t_eval=np.linspace(0, t_final, 100))
    return sol