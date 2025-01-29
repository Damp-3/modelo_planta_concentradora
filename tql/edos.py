##########################################
## PARAMETROS NO VARIABLES
##########################################

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from tql.parametros_variables import *
from scipy.integrate import quad

def pasteurizador_ode(t,y,params: Pastparams):
    """
    EDO para la temperatura del pasteurizador, Tpast(t).
    y = [Tpast]

    dTpast/dt = [ m_suero*cp_suero*(T_in - Tpast)
                  + U*A*(T_vapor - Tpast) ] / (m_hold * cp_suero)
    """
    Tpast = y[0]

    # extrae valores de params

    m_suero = params.m_suero
    f1 = params.f1
    T_in = params.get_T_in
    T_vap = params.get_T_vap
    cp_suero = params.get_cp_suero(T_in,f1)
    U = params.U
    A = params.get_area_past
    m_hold = params.get_m_hold(T_in,f1)


    Q_suero = m_suero * cp_suero * (T_in - Tpast)
    Q_vapor = U * A * (T_vapor - Tpast)
    
    dTpast_dt = (Q_suero + Q_vapor) / (m_hold * cp_suero)
    return [dTpast_dt]


