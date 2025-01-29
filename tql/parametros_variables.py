### parámetros que varian según la temperatura y/o el flujo másico

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from tql.propiedades import calcular_propiedades_agua
from tql.parametros_novariables import cargar_parametros_novariables
from scipy.integrate import quad

def h_cond(Tsat, T_w, param):
    """
    Calcula el coeficiente de transferencia de calor por condensación.
    
    Args:
        Tsat (float): Temperatura de saturación [K].
        T_w (float): Temperatura de la pared [K].
        param (dict): Parámetros no variables (L, nT, do).
        
    Returns:
        float: Coeficiente de transferencia de calor por condensación [W/m²·K].
    """
    # Obtener propiedades termodinámicas
    props = calcular_propiedades_agua(Tsat, T_w)
    
    # Propiedades necesarias
    k_l = props["conductividad_liquido"]
    rho_l = props["densidad_liquido"]
    visc_l = props["viscosidad_liquido"]
    rho_v = props["densidad_vapor"]
    h_lv_corr = props["entalpia_vapor_liquido"]
    
    # Parámetros geométricos y constantes
    L = param['L']
    nT = param['nT']
    do = param['do']
    g = param['g']  # Aceleración gravitacional [m/s²]
    b = math.pi * do  # Perímetro [m]
    
    # Fórmula del coeficiente de transferencia de calor
    h_cond =(0.943 * k_l / L) * ((rho_l * (rho_l - rho_v) * g * h_lv_corr * L**3) / (visc_l * k_l * (Tsat - T_w)))**(1 / 4)
    
    return h_cond


def h_conv(Tsat,T_w, mvap,efecto:int,param): 
    """
    Calcula el coeficiente de transferencia de calor por convección en tubos.
    
    Args:
        Tsat (float): Temperatura de saturación [K].
        mvap (float): Flujo de masa de vapor [kg/h].
        
    Returns:
        float: Coeficiente de transferencia de calor [W/m²·K].
    """
    # Obtener propiedades termodinámicas
    props = calcular_propiedades_agua(Tsat,T_w)
    
    # Propiedades específicas
    visc_v = props["viscosidad_vapor"]
    k_v = props["conductividad_vapor"]
    Prl = props["numero_prandtl_vapor"]
    
    if efecto == 1:
        d_evap = param['diametros_evaporadores']['1er efecto']
    elif efecto == 2:
        d_evap = param['diametros_evaporadores']['2do efecto']
    elif efecto == 3:
        d_evap = param['diametros_evaporadores']['3er efecto']
    elif efecto == 4:
        d_evap = param['diametros_evaporadores']['4to efecto']
    # Geometría de la tubería
    L = param["L"]  # Longitud [m]
    b = math.pi * d_evap  # Perímetro
    
    # Conversión de flujo de masa a kg/s
    mvap1 = mvap / 3600  # [kg/s]
    
    # Número de Reynolds
    Re = (4 * mvap1) / (visc_v * b)
    
    # Número de Nusselt (fórmula corregida)
    Nu = 0.0214 * (Re ** 0.8 - 100) * (Prl ** 0.4)
    
    # Coeficiente de transferencia de calor
    h_conv = Nu * k_v / L
    return h_conv


def h_evap(Tsat, T_w, m_2_e_p, param,efecto:int):
    """
    Calcula el coeficiente de transferencia de calor por evaporación.
    
    Args:
        Tsat (float): Temperatura de saturación [K].
        T_w (float): Temperatura de la pared [K].
        m_2_e_p (float): Flujo másico de evaporación por efecto [kg/h].
        param (dict): Parámetros no variables (L, di, nT).
        
    Returns:
        float: Coeficiente de transferencia de calor por evaporación [W/m²·K].
    """
    # Obtener propiedades termodinámicas
    props = calcular_propiedades_agua(Tsat, T_w)
    
    # Propiedades específicas
    k_l = props["conductividad_liquido"]
    rho_l = props["densidad_liquido"]
    visc_l = props["viscosidad_liquido"]
    
    # Parámetros geométricos
    L = param["L"]
    di = param["di"]
    nT = param["nT"]
    g = param["g"]  # Aceleración gravitacional [m/s²]
    
    # Flujo másico unitario por tubo [kg/s]
    if efecto == 1:
        nT = param['tubos_por_efecto']['1er efecto']
    elif efecto == 2:
        nT = param['tubos_por_efecto']['2do efecto']
    elif efecto == 3:
        nT = param['tubos_por_efecto']['3er efecto']
    elif efecto == 4:
        nT = param['tubos_por_efecto']['4to efecto']

    m_2_e_p_unit = m_2_e_p / (3600 * nT)
    
    # Número de Reynolds
    Re = (4 * m_2_e_p_unit) / (visc_l * math.pi * di)
    
    # Cálculo del coeficiente de transferencia de calor por evaporación
    h_evap = (((4 * rho_l * g * k_l**3) / (3 * visc_l**2))**(1/3)) * Re**(-1/3)
    
    return h_evap


def rho_prod(T,f1):
    """Calcula la densidad del producto en función de la temperatura y el porcentaje de solidos totales.
    
    Parameters
    ----------
    T : float
        Temperatura del producto en °C
    f1 : float
        Porcentaje de solidos totales
    """
    T_k = T+273
    rho_agua = PropsSI('D','T',T_k,'Q',0,'water')
    rho_prod = (1100 - rho_agua *(1-f1))/f1
    return rho_prod

def visc_suero(T,f1):
    """Calcula la viscosidad del suero de la leche en función de la temperatura y el contenido de sólidos totales
    Parameters
    ----------
    T : float
        Temperatura del producto en °C
    f1 : float   
        Contenido en porcentaje de solidos totales valor entre 0 y 1
    Returna:
    - viscosidad(float): viscosidad del suero de la leche en Pa.s
    """
    # Parámetros ajustables del modelos
    A = 1.1 # mPa.s a 0°C
    B = 0.045 # Factor de disminión con la temperatura
    C = 1.5 # Ajuste según contenido de agua

    # Efecto de la temperatura 
    factor_temp = np.exp(-B*T)
    # Efecto del contenido de agua
    factor_agua = f1**C
    # Viscosidad resultante
    viscosidad = A*factor_temp*factor_agua
    return viscosidad*1000



def cp_suero(T, f1):
    """
    Calcula la capacidad calorífica del suero de leche en función de:
      - contenido de sólidos totales (f1, fracción entre 0 y 1),
      - temperatura (T, en °C).
    Retorna la cp en J/kg.K (o J/kg.°C, equivalentes para incrementos de T).
    """

    # Asegurarnos de convertir T a Kelvin si usamos CoolProp:
    T_K = T + 273.15

    # Ejemplo: cp del agua usando CoolProp en J/kg.K a la temperatura dada:
    cp_agua_0 = PropsSI('C', 'T', T_K, 'Q', 0, 'Water')  # entalpía mas. en J/kg.K (para agua líquida)
    
    # Valores para sólidos totales, etc. (según tu modelo, aquí en un ejemplo simple)
    cp_solid_0 = 1800    # J/kg.K (sólidos)
    
    # Coeficientes de ajuste por temperatura
    k_agua = 0.001       # ejemplo de incremento por °C para agua
    k_solid = 0.0005     # ejemplo de incremento por °C para sólidos
    
    # Ajustamos con T (en °C):
    cp_agua = cp_agua_0 + k_agua * T
    cp_solid = cp_solid_0 + k_solid * T

    # Fracción de agua
    contenido_agua = 1.0 - f1

    # Mezcla ponderada
    cp = cp_agua * contenido_agua + cp_solid * f1
    return cp

def Q_suero_integral(m_suero, T_in, T_out, f1):
    """
    Calcula la energía necesaria (en J/h) para calentar el flujo de suero 
    desde T_in hasta T_out, integrando cp_suero(T, f1).
    
    m_suero_h : float
        Flujo másico de suero en kg/h.
    T_in : float
        Temperatura inicial (°C).
    T_out : float
        Temperatura final (°C).
    f1 : float
        Fracción de sólidos totales (0-1).
    """
    # Función que integrará cp_suero(T, f1) con respecto a T
    def integrand(T):
        return cp_suero(T, f1)  # [J/kg.K]

    # Integración de cp_suero entre T_in y T_out
    # Esto nos dará algo en [J/kg] (pues integras J/kg.K * dT)
    Q,_ = quad(integrand, T_in, T_out)

    # Ahora multiplicamos por el caudal másico en kg/h, y dividimos por 3600 para que de en [W]
    Q_total = m_suero * Q /3600

    return Q_total

def Delta_h_subenf(Tsat,T):
    """Calcula la entalpía de vaporización del agua subenfriado a temperatura T < Tsat

    Parameters
    ----------
    Tsat : float
        Temperatura de saturación °C
    T : float
        Temperatura de subenfriamiento °C
    Retorna: Entalpía de vaporización del agua a temperatura T en [J/kg]
    """
    Tsat_k = Tsat + 273
    T_k = T + 273
    cp_l = PropsSI('C', 'T', Tsat_k, 'Q',0,'water')
    h_l = PropsSI('H', 'T', Tsat_k, 'Q',0,'water')
    h_v = PropsSI('H', 'T', Tsat_k, 'Q',1,'water')
    h_lv = h_v-h_l

    h_lv_sub = h_lv + (Tsat_k-T_k)*cp_l

    return h_lv_sub

def h_l(T,f1):
    """Calcula la entalpía del suero de leche en función de la temperatura y los sólidos totales.

    Parameters
    ----------
    T : float
        Temperatura actual del suero de leche en [°C]
    f1 : float
        Contenido de sólidos totales
    tref : float    
        Temperatura de referencia en [°C]
    """

    h_l, _ = quad(lambda T: cp_suero(T,f1),0,T)
    return h_l  

#m_suero,m_hold,T_in,T_vap,A
class Pastparams:
    def __init__(self):
        self.m_suero = np.random.uniform(19000,21000)
        self.L = 13
        self.e = 0.005
        self.g = 9.81
        self.do = 0.0508
        self.num_tub = 2
        

        def area_past(self):
            return math.pi*self.L*self.do*self.num_tub

        def di(self):
            return self.do - 2 * self.e
        def rho_prod(self,T,f1):
            """Calcula la densidad del producto en función de la temperatura y el porcentaje de solidos totales.
            Parameters
            ----------
            T : float
                Temperatura del producto
            f1 : float
                Porcentaje de solidos totales
            """
            T_k = T+273
            rho_agua = PropsSI('C','T',T_k,'Q',1,'water')
            rho_prod = (1.10-rho_agua*(1-f1))/f1
            return rho_prod



        def m_hold(self):
            return math.pi
        pass
