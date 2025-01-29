### parámetros que varian según la temperatura y/o el flujo másico

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from tql.propiedades import calcular_propiedades_agua
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


class Pastparams:
    """
    Clase para encapsular parámetros y métodos físicos del pasteurizador.
    
    Atributos (pasados al constructor):
    ----------------------------------
    L        : Longitud total de los tubos (m)
    e        : Espesor de los tubos (m)
    g        : Aceleración gravitatoria (m/s^2)
    do       : Diámetro externo del tubo (m)
    num_tub  : Número de tubos
    f1       : Fracción (0 a 1) de sólidos totales
    """
    def __init__(self, 
                 L: float = 13.0, 
                 e: float = 0.005, 
                 g: float = 9.81, 
                 do: float = 0.0508, 
                 num_tub: int = 2, 
                 f1: float = 0.115,
                 U: float = 7000):
        
        self.L = L
        self.e = e
        self.g = g
        self.do = do
        self.num_tub = num_tub
        self.f1 = f1
        self.U = U
    
    def m_suero(self):
        """
        Retorna un flujo de suero aleatorio entre 19000 y 21000 (unidades a definir).
        
        """
        return np.random.uniform(19000, 21000)

    def get_area_past(self):
        """
        Retorna el área externa total de los tubos,
        asumiendo do como diámetro externo, L como longitud del tubo,
        y num_tub como cantidad de tubos.
        A = pi * do * L * num_tub
        """
        return math.pi * self.L * self.do * self.num_tub

    def get_di(self):
        """
        Diámetro interno: do - 2*e
        """
        return self.do - 2.0*self.e
    
    def get_rho_prod(self, T: float, f1: float):
        """
        Calcula la densidad del producto en función de la temperatura (°C)
        y la fracción de sólidos totales (f1).
        
        Usamos CoolProp para la densidad del agua a T dada (en K),
        y sumamos la contribución de sólidos (densidad fija 1397 kg/m3).
        """
        T_k = T + 273.15
        rho_agua = PropsSI('D','T', T_k,'Q',0,'Water')  # densidad del agua líquida
        rho_prod = rho_agua*(1 - f1) + 1397*f1
        return rho_prod

    def get_m_hold(self, T: float, f1: float):
        """
        Retorna la masa de producto en el interior (hold-up),
        asumiendo un volumen cilíndrico (sección interna) * longitud,
        multiplicado por la densidad del producto (get_rho_prod).
        """
        di = self.get_di()
        area_interna = math.pi * (di**2) / 4.0  # sección transversal
        # Ojo: en tu código pusiste math.pi * (di**2), que es el área de un círculo con radio=di. 
        # Normalmente, area = pi*(diametro^2)/4. Ajusta si lo tuyo es otra geometría.
        
        rho = self.get_rho_prod(T, f1)
        volume = area_interna * self.L
        return volume * rho

    def get_cp_suero(self, T: float, f1: float):
        """
        Calcula la capacidad calorífica del suero (J/kg.K), 
        suponiendo una mezcla de agua + sólidos.
        
        - T (°C) => se convierte a Kelvin para usar CoolProp.
        - f1 = fracción de sólidos (0..1).
        
        Se obtiene cp_agua_0 de CoolProp y se corrige levemente con un factor.
        Se usa un cp_solid_0 aproximado para la parte sólida.
        """
        T_K = T + 273.15
        # cp_agua_0 => calor específico del agua en J/kg.K
        cp_agua_0 = PropsSI('C', 'T', T_K, 'Q', 0, 'Water')
        
        # supuestos para los sólidos
        cp_solid_0 = 1800.0   # J/kg.K (sólidos)
        
        # "ajustes" por temperatura (ejemplo hipotético)
        k_agua = 0.001   # J/kg.K por °C
        k_solid = 0.0005
        
        cp_agua = cp_agua_0 + k_agua * T
        cp_solid = cp_solid_0 + k_solid * T
        
        # Mezcla ponderada
        contenido_agua = 1.0 - f1
        cp_mixture = cp_agua*contenido_agua + cp_solid*f1
        
        return cp_mixture

    def get_T_in(self):
        """Genera una temperatura aleatoria de ingreso de producto a partir de una distribución uniforme entre 70 y 74 °C

        Returns
        -------
        float
            Temperatura de ingreso de producto en °C
        """
        return np.random.uniforme(70,74)
    
    def get_T_vap(self):
        """Genera una temperatura aleatoria de ingreso de vapor a partir de una distribución uniforme entre 77 y 83 °C

        Returns
        -------
        float
            Temperatura de ingreso de producto en °C
        """
        return np.random.uniforme(77,83)
    
