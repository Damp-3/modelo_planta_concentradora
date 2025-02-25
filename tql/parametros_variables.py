### parámetros que varian según la temperatura y/o el flujo másico

import math
import numpy as np
from CoolProp.CoolProp import PropsSI
from scipy.integrate import quad


def visc_suero(T,f1):
    """Calcula la viscosidad del suero de la leche en función de la temperatura y el contenido de sólidos totales
    Parameters
    ----------
    T : float
        Temperatura del producto en K
    f1 : float   
        Contenido en porcentaje de solidos totales valor entre 0 y 1
    Returna:
    - viscosidad(float): viscosidad del suero de la leche en Pa.s
    """
    visc_agua = PropsSI('V','T',T,'Q',0,'water')
    m = 1
    n = 3.5
    K = 8
    visc_s = visc_agua*(1 + (K*f1)**n)**m
    return visc_s

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
    U        : Coeficiente global de transferencia de calor (W/m2/K)
    """
    def __init__(self, 
                 L: float = 13.0, 
                 e: float = 0.005, 
                 g: float = 9.81, 
                 do: float = 0.0508, 
                 num_tub: int = 2, 
                 f1: float = 0.115,
                 U: float = 4000):
        
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
    
    def get_rho_prod(self, T_K: float, f1: float):
        """
        Calcula la densidad del producto en función de la temperatura (°C)
        y la fracción de sólidos totales (f1).
        
        Usamos CoolProp para la densidad del agua a T dada (en K),
        y sumamos la contribución de sólidos (densidad fija 1397 kg/m3).
        """
        rho_agua = PropsSI('D','T', T_K,'Q',0,'Water')  # densidad del agua líquida
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

    def get_cp_suero(self, T_K: float, f1: float):
        """
        Calcula la capacidad calorífica del suero (J/kg.K), 
        suponiendo una mezcla de agua + sólidos.
        
        - T (K) => se convierte a Kelvin para usar CoolProp.
        - f1 = fracción de sólidos (0..1).
        
        Se obtiene cp_agua_0 de CoolProp y se corrige levemente con un factor.
        Se usa un cp_solid_0 aproximado para la parte sólida.
        """
        # cp_agua_0 => calor específico del agua en J/kg.K
        cp_agua_0 = PropsSI('C', 'T', T_K, 'Q', 0, 'Water')
        
        # supuestos para los sólidos
        cp_solid_0 = 1800.0   # J/kg.K (sólidos)
        
        # "ajustes" por temperatura (ejemplo hipotético)
        k_agua = 0.001   # J/kg.K por °C
        k_solid = 0.0005
        
        cp_agua = cp_agua_0 + k_agua * T_K
        cp_solid = cp_solid_0 + k_solid * T_K
        
        # Mezcla ponderada
        contenido_agua = 1.0 - f1
        cp_mixture = cp_agua*contenido_agua + cp_solid*f1
        
        return cp_mixture

    def get_T_in(self):
        """Genera una temperatura aleatoria de ingreso de producto a partir de una distribución uniforme entre 70 y 74 °C

        Returns
        -------
        float
            Temperatura de ingreso de producto en K
        """
        return np.random.uniform(343,347)
    
    def get_T_vap(self):
        """Genera una temperatura aleatoria de ingreso de vapor a partir de una distribución uniforme entre 77 y 83 °C

        Returns
        -------
        float
            Temperatura de ingreso de vapor en K
        """
        return np.random.uniform(350,356)
    
class EvapEffectParams:
    """
    Clase para modelar los efectos de evaporación en un sistema de evaporación múltiple.

    Atributos:
    ----------
    Pastparams : object
        Objeto que contiene parámetros generales del sistema.
    efecto : int
        Número del efecto en el sistema de evaporación.
    L : float
        Longitud del tubo [m].
    g : float
        Aceleración gravitatoria [m/s²].
    do : float
        Diámetro externo del tubo [m].
    f1 : float
        Fracción de sólidos totales.
    num_tub : int
        Número de tubos en el efecto.
    d_evap : float
        Diámetro del evaporador en el efecto.
    """
    def __init__(self, Pastparams, efecto: int):
        self.Pastparams = Pastparams  # Parámetros generales del sistema
        self.efecto = efecto          # Número del efecto
        self.L = 13.0                 # Longitud del tubo [m]
        self.g = 9.81                 # Aceleración gravitatoria [m/s²]
        self.do = 0.0508              # Diámetro externo del tubo [m]
        self.di = 0.0488              # Diámetro interno del tubo [m]
        self.f = {1: 0.115, 2: 0.158, 3: 0.25, 4: 0.332}[efecto] # Fracción de sólidos totales por efecto
        self.num_tub = {1: 180, 2: 180, 3: 120, 4: 120}[efecto]  # Tubos por efecto
        self.d_evap = {1: 1.00, 2: 1.00, 3: 0.95, 4: 0.95}[efecto]  # Diámetro del evaporador
        self.cp_inox = 500                                          # Capacidad calorífica del inoxidable [J/kgK]
        self.d_inox = 7800                                          # Densidad acero inoxidable [kg/m3]
    

    def m_inox(self):
        """Define la masa total de un tubo
        """
        m_inox = (math.pi*((self.do/2)**2 - (self.di/2)**2))*self.L * self.d_inox

        return m_inox

    def get_area_effect(self,locus:str):
        """Retorna el área de intercambio interna o externa de 1 tubo en [m2].

        Parameters
        ----------
        locus : str
            "in" retorna el área interna o "out" retorna el área externa
        efecto : int
            n° de efecto relacionado
        """

        if locus == "in":
            A = math.pi*self.di * self.L
        else:
            A = math.pi*self.do * self.L

        return A
    


    def ent_latent_cond(self, Ts, Tw):

        Tw = np.array(Tw) # Convertir a array si no lo es
        T_prom = (Ts + Tw) / 2

        # Evaluar PropsSI solo en valores escalares
        cp_agua = PropsSI('C','T',Ts,'Q',0,'water')

        h_l = PropsSI('H', 'T', Ts, 'Q', 0, 'water')
        h_v = PropsSI('H', 'T', Ts, 'Q', 1, 'water')
        h_lv_corr = (h_v - h_l) + (0.68 * cp_agua * (Ts - Tw))

        # Validar h_lv_corr
        if np.isnan(h_lv_corr) or h_lv_corr <= 0:
            print(f"Advertencia: h_lv_corr inválido ({h_lv_corr}). Ajustando...")
            h_lv_corr = 2.328e6  # Valor seguro

        return h_lv_corr

    def h_cond(self, Tsat, T_w):
        """
        Calcula el coeficiente de transferencia de calor por condensación.

        Parámetros:
        -----------
        Tsat : float
            Temperatura de saturación [K].
        T_w : float
            Temperatura de la pared [K].

        Retorna:
        --------
        float
            Coeficiente de transferencia de calor por condensación [W/m²·K].
        """

        if T_w >= Tsat:
        # La pared está más caliente que el vapor => no hay condensación
            return 0
        Temp_pel = (Tsat + T_w)/2
        k_l = PropsSI('L', 'T', Temp_pel, 'Q', 0, 'water')             # [W/m·K]
        rho_l = PropsSI('D', 'T', Temp_pel, 'Q', 0, 'water')           # [kg/m³]
        visc_l = PropsSI('V', 'T', Temp_pel, 'Q', 0, 'water')          # [Pa.s]
        rho_v = PropsSI('D', 'T', Temp_pel, 'Q', 1, 'water')           # [kg/m³]
        h_lv_corr = self.ent_latent_cond(Tsat,T_w)
        Prl = PropsSI('Prandtl','T',Temp_pel,'Q',0,'water')
        g = self.g
        L = self.L

        # Validar delta_T
        delta_T = Tsat - T_w
        if delta_T <= 0:
            print(f"Advertencia: Tsat ({Tsat}) - T_w ({T_w}) no es válido. Ajustando...")
            delta_T = max(delta_T, 1e-3)

        # Validar h_lv_corr antes de usarlo
        if np.isnan(h_lv_corr) or h_lv_corr <= 0:
            print(f"Advertencia: h_lv_corr inválido ({h_lv_corr}). Ajustando...")
            h_lv_corr = 2.3e6  # Valor seguro

        # Calcular h_cond evitando valores negativos o NaN
        try:
            h_cond = 0.943 * ((rho_l * (rho_l - rho_v) * g * h_lv_corr * k_l**3) / (visc_l * delta_T * L))**(1/4)
        except ValueError:
            print("Error en el cálculo de h_cond, ajustando valores...")
            h_cond = 2200  # Valor seguro

        if np.isnan(h_cond) or h_cond < 0:
            print(f"Advertencia: h_cond inválido ({h_cond}). Ajustando a 2200.")
            h_cond = 2200  # Valor seguro

        return h_cond
    
    def h_conv(self, Tsat, T_w, mvap):
        """
        Calcula el coeficiente de transferencia de calor por convección en tubos.

        Parámetros:
        -----------
        Tsat : float
            Temperatura de saturación [°C].
        T_w : float
            Temperatura de la pared [°C].
        mvap : float
            Flujo de masa de vapor [kg/h].

        Retorna:
        --------
        float
            Coeficiente de transferencia de calor por convección [W/m²·K].
        """
        visc_v = PropsSI('V','T',Tsat,'Q',1,'water')
        k_v = PropsSI('L','T',Tsat,'Q',1,'water')
        Prl = PropsSI('Prandtl','T',Tsat,'Q',1,'water')
        
        b = math.pi * self.d_evap  # Perímetro del tubo
        mvap1 = mvap / 3600  # Conversión de flujo de masa a kg/s
        
        Re = (4 * mvap1) / (visc_v * b)
        Nu = 0.0214 * (Re ** 0.8) * (Prl ** 0.4)
        
        h_conv = Nu * k_v / self.L
        return h_conv
    
    def h_evap(self, Tsat, T_w, m_2_e_p):
        """
        Calcula el coeficiente de transferencia de calor por evaporación.

        Parámetros:
        -----------
        Tsat : float
            Temperatura de saturación [K].
        T_w : float
            Temperatura de la pared [K].
        m_2_e_p : float
            Flujo másico de evaporación por efecto [kg/h].

        Retorna:
        --------
        float
            Coeficiente de transferencia de calor por evaporación [W/m²·K].
        """
        T_prod = self.Pastparams.get_T_in()
        f = self.f
        k_l = PropsSI('L','T',Tsat,'Q',0,'water')
        visc_s = visc_suero(T_prod,f)
        rho_s = self.Pastparams.get_rho_prod(T_prod,f)
        g = self.g


        m_2_e_p_unit = m_2_e_p / (3600 * self.num_tub)
        
        Re = (4 * m_2_e_p_unit) / (visc_s * math.pi * self.di)
        
        h_evap = (((4 * rho_s**2 * g * k_l**3) / (3 * visc_s**2))**(1/3)) * Re**(-1/3)
        return h_evap
    
    def m_vap_1(self,x):
        """Retorna según efecto la masa de vapor que pasa por el volumen de control de cada efecto.
        Se hacen alguna suposiciones de cuánto vapor condensa en cada efecto.
        Parameters
        ----------
        efecto : int
            Efecto en el cual se está.
        Retorna:
            El valor (float) de el flujo másico del vapor en [kg/h].
        """ 
        # se asume a la salida que sólo un 85% de flujo de vapor sigue siendo efectivamente vapor saturaado.
        m_vap_in = x # tasa de vapor de entrada efecto 1 a 80°C desde manifold en [kg/h]
        m_vap_out = m_vap_in*0.85 # tasa de vapor de salida efecto 1 a 80°C 
        return m_vap_in,m_vap_out
    
















        # elif efecto == 2:
        #     # se asume a la salida que sólo un 85% de flujo de vapor sigue siendo efectivamente vapor saturaado.
        #     m_evap_gen = self.m_evap(1) # tasa de vapor generado por evaporación (SS) [kg/h] 
        #     m_vap_in = 1400 + m_evap_gen # flujo de vapor que entra al efecto 2.
        #     m_vap_out = (m_vap_in)*0.85 # tasa de vapor que sale del efecto 1, y va a la entrada del efecto 2.
        #     return m_vap_in, m_vap_out
        # elif efecto == 3:
        #     # se asume a la salida que sólo un 85% de flujo de vapor sigue siendo efectivamente vapor saturaado.
        #     m_evap_gen = self.m_evap(2) # tasa de vapor generado por evaporación (SS) [kg/h]
        #     m_vap_in = ((1400 + self.m_evap(1))*0.85 + self.m_evap(2))*0.85 # tasa de vapor que entra al efecto 3
        #     m_thermo_c = m_vap_in*0.5*0.85
        #     m_vap_out = ((1400 + self.m_evap(1))*0.85 + self.m_evap(2))*0.85 + self.m_evap(2)*0.5 * 0.85 # tasa de vapor que sale del efecto 2, y va a la entrada del efecto 3.
        #     return m_vap_in, m_vap_out
        # else:
        #     # se asume a la salida que sólo un 85% de flujo de vapor sigue siendo efectivamente vapor saturaado.
        #     m_evap_gen = self.m_evap(3) # tasa de vapor generado por evaporación (SS) [kg/h]
        #     m_vap_in = (1400 + ((1400 + self.m_evap(1))*0.85 + self.m_evap(2))*0.85 + self.m_evap(2)) * 0.85 # tasa de vapor que entra al efecto 3
        #     m_vap_out =

def ciclo_m_vap():
    results = [[],[],[],[]]
    for i in range(1,4):

        m_i, m_o = m_vap(i)
        r = [m_i, m_o]
        results[i] = r

    return results


