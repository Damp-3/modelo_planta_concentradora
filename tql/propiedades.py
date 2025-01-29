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

def calcular_propiedades_agua(Tsat, T_w=None ):
    """
    Calcula propiedades termodinámicas necesarias para h_cond y devuelve un diccionario.
    
    Args:
        Tsat (float): Temperatura de saturación [K].
        T_w (float): Temperatura de la pared [K] (opcional).
        
    Returns:
        dict: Diccionario con las propiedades calculadas.
    """
    h_l = PropsSI('H', 'T', Tsat, 'Q', 0, 'water')  # [J/kg]
    h_v = PropsSI('H', 'T', Tsat, 'Q', 1, 'water')  # [J/kg]
    cp_l = PropsSI('C', 'T', Tsat, 'Q', 0, 'water')  # [J/kg·K]
    h_lv = h_v - h_l  # [J/kg]
    
    # Corrección de entalpía si T_w está disponible
    h_lv_corr = h_lv + 0.68 * cp_l * (Tsat - T_w) if T_w else h_lv
    
    return {
        "conductividad_liquido": PropsSI('L', 'T', Tsat, 'Q', 0, 'water'),  # [W/m·K]
        "densidad_liquido": PropsSI('D', 'T', Tsat, 'Q', 0, 'water'),       # [kg/m³]
        "viscosidad_liquido": PropsSI('V', 'T', Tsat, 'Q', 0, 'water'),     # [Pa.s]
        "densidad_vapor": PropsSI('D', 'T', Tsat, 'Q', 1, 'water'),         # [kg/m³]
        "capacidad_calorifica_liquido": cp_l,                              # [J/kg·K]
        "entalpia_liquido": h_l,                                           # [J/kg]
        "entalpia_vapor": h_v,                                             # [J/kg]
        "numero_prandtl_liquido": PropsSI('Prandtl', 'T', Tsat, 'Q', 0, 'water'),  # [-]
        "entalpia_vapor_liquido": h_lv_corr,                               # [J/kg]   
        "conductividad_vapor":PropsSI('L','T', Tsat,'Q',1,'water'),
        "viscosidad_vapor": PropsSI('V','T',Tsat,'Q',1,'water'),
        "capacidad_calorifica_vapor": PropsSI('C','T',Tsat,'Q',1,'water'),
        "numero_prandtl_vapor": PropsSI('Prandtl', 'T', Tsat, 'Q', 1, 'water')   # [J/kg.K]

    }

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