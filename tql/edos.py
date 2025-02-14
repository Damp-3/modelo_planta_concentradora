##########################################
## PARAMETROS NO VARIABLES
##########################################

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from tql.parametros_variables import *
from scipy.integrate import quad
from scipy.integrate import solve_ivp


def pasteurizador_ode(t,y,params: Pastparams):
    """
    EDO para la temperatura del pasteurizador, Tpast(t).
    y = [Tpast]

    dTpast/dt = [ m_suero*cp_suero*(T_in - Tpast)
                  + U*A*(T_vapor - Tpast) ] / (m_hold * cp_suero)
    """
    Tpast = y[0]

    # extrae valores de params

    m_suero = params.m_suero()
    f1 = params.f1
    T_in = params.get_T_in()
    T_vap = params.get_T_vap()
    cp_suero = params.get_cp_suero(T_in,f1)
    U = params.U
    A = params.get_area_past()
    m_hold = params.get_m_hold(T_in,f1)


    Q_suero = m_suero * cp_suero * (T_in - Tpast)
    Q_vapor = U * A * (T_vap - Tpast)
    
    dTpast_dt = (Q_suero + Q_vapor) / (m_hold * cp_suero)
    return [dTpast_dt]

def simulate_pasteurization(t_span, T0, params: Pastparams):
    """
    Integra la ODE del pasteurizador en [t_span].
    T0: temperatura inicial.
    params: instancia de PasteurParams.
    """
    sol = solve_ivp(
        fun=lambda t, y: pasteurizador_ode(t, y, params),
        t_span=t_span,
        y0=[T0],
        dense_output=True
    )
    return 

def evap_ode(t,T, pastparams: Pastparams,params:EvapEffectParams):
    """Define la ecuación diferencial para la evolución de la temperatura de la pared
    en el primer efecto de evaporación.

    dTdt = (h_cond * A * (T_vapor - T)-h_evap * A * (T - T_product))/(m_wall * c_wall)
    """
    effectparams_e1 = EvapEffectParams(None,1)
    Tvap = Pastparams.get_T_vap()
    Tprod = 73                       #sol.simulate_pasteur
    m_prod_in = Pastparams.m_suero()
    h_cond = EvapEffectParams.h_cond(Tvap,T)
    h_evap = EvapEffectParams.h_evap(Tvap,T,m_prod_in)
    A_in = EvapEffectParams.get_area_effect("in")
    A_out = EvapEffectParams.get_area_effect("out")
    m_inox = EvapEffectParams.m_inox()
    cp_inox = EvapEffectParams.cp_inox()

    Q_cond = h_cond * A_out * (Tvap - T)
    Q_evap  = h_evap * A_in * (T - Tprod)

    dTdt = (Q_cond - Q_evap) / (m_inox * cp_inox)
    return dTdt

def sim_evap(T0,t_final,params,evap_effect):
    """Simula la ecuación diferencial para el primer efecto de evaporación.

    Parameters
    ----------
    T0 : float
        Temperaturas iniciales en [°C]
    t_final : float
        Tiempo de simulación en [s]
    params : clase
        Clase con atributos y métodos de la pasteruización
    evap_effect : clase
        Clase con atributos y métodos del efecto.
    """

    sol_e1 = solve_ivp(lambda t, T: evap_ode(t, T, params, evap_effect),
                       [0, t_final], [T0], t_eval = np.linspace(0,t_final,100))
    
    return sol_e1
