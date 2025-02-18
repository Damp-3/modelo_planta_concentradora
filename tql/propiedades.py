import pandas as pd
from CoolProp.CoolProp import PropsSI
import numpy as np
import math
from tqdm import tqdm
from CoolProp.CoolProp import AbstractState
import CoolProp

def h_lv_corr(Tsat1,T):
    h_l = PropsSI('H', 'T', Tsat1, 'Q', 0, 'water') # [J/kg]
    h_v = PropsSI('H', 'T', Tsat1, 'Q', 1, 'water') # [J/kg]
    h_lv = h_v - h_l # [J/kg]
    cp_l = PropsSI('C', 'T', Tsat1, 'Q', 0, 'water') # [J/kg/K]
    h_lv_c = h_lv + 0.68*cp_l*(T - Tsat1) # [J/kg]
    return h_lv_c

def cp_air(T,P):
    """Calcula la capacidad calorífica del aire a diferentes Temperaturas

    Parameters
    ----------
    T : float
        Temperatura del aire en °C
    """
    
    HEOS = CoolProp.AbstractState("HEOS", "Air")
    HEOS.update(CoolProp.PT_INPUTS, P, T)
    cp_air = HEOS.cpmass()
    return cp_air