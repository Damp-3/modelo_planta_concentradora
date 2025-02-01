######################################################
## Modelo 1er efecto
######################################################
from CoolProp.CoolProp import PropsSI
import numpy as np
import pandas as pd
from tql.parametros_variables import *  # noqa: F403
from tql.propiedades import cp_air
from tql.utils import *  # noqa: F403
from scipy.integrate import quad
from tql.modelos import *  # noqa: F403
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


def entre_efecto(Tsat1_E1,Tsat2_E1,T_w):
    """Esta función relaciona la salida de un efecto, con el inicio del siguiente efecto
    Returns
    -------
    tuple
        Retorna, las temperaturas y presiones que debe tener el inicio del siguiente efecto.
    """
    Tsat1_E2 = (Tsat1_E1+T_w)/(2)
    Tsat2_E2 = (Tsat2_E1+T_w)/(2)
    P1_E2 = PropsSI('P','T',Tsat1_E2,'Q',0.5,'water')
    P2_E2 = PropsSI('P','T',Tsat2_E2,'Q',0.5,'water')



    return Tsat1_E2, Tsat2_E2, P1_E2, P2_E2

# fors



def efecto_estacionario(efectos):
    """Esta función calcula el flujo total de vapor requerido para evaporar en estado estacionario.
    Parameters
    ----------
    efectos : N° de efectos
    """
    # _,b = cargar_datos()
    # m_1_e_v = a
    # m_2_e_p = b
    Tsat1 = np.array([80,70,60,50])+273
    T_intermedia = np.array([70,63,55,48])+273
    Tsat2 = np.zeros(len(Tsat1))
    f = [0.115,0.25,0.332,0.48]

    for i in range(len(Tsat1)):
        Tsat2[i]= (T_intermedia[i]+Tsat1[i])/2

    m_evap = [6160,5930,2597,2433]
    m_vap_req = np.zeros(len(Tsat1))


    for efecto in range(efectos):
        # lado producto
        h_l_p_E1 = h_l(Tsat2[efecto],f[efecto])
        h_v_p_E1 = PropsSI('H','T',Tsat2[efecto],'Q',1,'water')
        h_lv_prod = h_v_p_E1 - h_l_p_E1 # en [J/kg]
        # lado vapor
        h_l_v_E1 = PropsSI('H','T',Tsat1[efecto],'Q',0,'water')
        h_v_v = PropsSI('H','T',Tsat1[efecto],'Q',1,'water')
        h_lv_vap = h_v_v - h_l_v_E1 # en [J/kg].

        m_vap_req[efecto] = m_evap[efecto]*h_lv_prod/h_lv_vap

    m_vap_req_t = sum(m_vap_req)
    return m_vap_req_t


def KPI_evap(efectos,m_evap_t):
    """Calcula el ratio entre kg de vapor/ kg de producto

    Returns
    -------
    float
        Valor del ratio entre kg de vapor /kg de producto procesado
    """
    m_vap_total = efecto_estacionario(efectos) 

    kpi = m_vap_total/m_evap_t
    return kpi



# def efecto_no_estacionario(Tsat1_E1,Tsat2_E2,T_w_E1,efecto:int):
#     """Esta función 

#     Parameters
#     ----------
#     Tsat1_E1 : _type_
#         _description_
#     Tsat2_E2 : _type_
#         _description_
#     T_w_E1 : _type_
#         _description_
#     """
#     T_w_E1 = 310.15 # 75°C
#     Tsat1_E1 = range(313,353,1) # 40-80°C
#     Tsat2_E1 = range(303,343,1) # 30-70°C
#     Tamb = 293.15 # 20°C
#     N = len(mv_evap)
#     M = len(Tsat1_E1)
#     L = len(Tsat2_E1)
#     # masas
#     m_1_e_v_E1 = mv_evap
#     m_1_s_v_E1 = np.zeros(N) # init
#     m_1_s_c_E1 = np.zeros(N) # init
#     m_e_p_E1 = np.mean(m_2_e_p)
#     m_2_s_p_E1 = np.zeros(N) # init
#     m_2_s_v_E1 = np.zeros(N) # init
#     #fluxos de calor
#     q_1_e_v_E1 = np.zeros((N,M)) # init
#     q_1_i_E1 = np.zeros((N,M)) # init
#     q_1_s_v_E1 = np.zeros((N,M)) # init
#     q_1_s_c_E1 = np.zeros((N,M)) # init
#     q_2_s_p_E1 = np.zeros((N,L)) # init
#     q_2_i_E1 = np.zeros((N,L)) # init
#     q_2_s_i_E1 = np.zeros((N,L)) # init
#     q_1_s_a_E1 = np.zeros(M)
#     q_1_s_w_E1 = np.zeros(M)
#     q_2_e_p_E1 = np.zeros(L)


#     for i in tqdm(range(N)):
        
#         m_1_s_v_E1[i] = m_1_e_v_E1[i] - m_1_i_E1[i]
#         m_1_s_c_E1[i] = m_1_i_E1[i]
#         m_2_s_p_E1[i] = m_2_e_p[i] - m_2_i_E1[i]
#         m_2_s_v_E1[i] = m_2_i_E1[i]

#         # Calores entrada y salida de vapor
#         for j in range(M):
#             q_1_e_v_E1[i,j] = m_1_e_v_E1[i] * cp_v(Tsat1_E1[j]) * Tsat1_E1[j]/3600
#             q_1_s_a_E1[j] = h_conv(Tsat1_E1[j],mv_evap[i]) * (Tsat1_E1[j] - Tamb)*a_l_evap_1
#             q_1_i_E1[i,j] = m_1_i_E1[i] * h_lv_corr(Tsat1_E1[j],T_w)/3600
#             q_1_s_v_E1[i,j] = m_1_s_v_E1[i] * PropsSI('C', 'T', Tsat1_E1[j],'Q', 1, 'water') * j/3600

#         # q_1_e_v[i] = q_1_s_a + q_1_i[i] + q_1_s_v[i]

#         # Calores entrada y salida de condensado

#             q_1_i_E1[i,j] = m_1_i_E1[i] * h_lv_corr(Tsat1_E1[j],T_w)/3600
#             q_1_s_w_E1[j] = h_cond(Tsat1_E1[j],T_w) * (Tsat1_E1[j] - T_w) * alt_1
#             q_1_s_c_E1[i,j] = m_1_s_c_E1[i] * cp_l(Tsat1_E1[j] )*  Tsat1_E1[j] / 3600
        
#             # q_1_i[i] = q_1_s_w + q_1_s_c[i]
        
#         # Calores entrada y salida producto (suero)
#             for l in range(L):
#                 q_2_e_p_E1[l] = h_evap(Tsat2_E1[l],T_w,m_e_p) * (Tsat2_E1[l] - T_w) * alt_1
#                 q_2_s_p_E1[i,l] = m_2_s_p_E1[i] * cp_l(Tsat2_E1[l]) * Tsat2_E1[l]/3600
#                 q_2_i_E1[i,l] = m_2_i_E1[i] * h_lv_corr(Tsat2_E1[l],T_w)/3600
            
#                 # q_2_e_p = q_2_s_p[i] + q_2_i[i]
            
#             # Calores entrada y salida de vapor producto
            
#                 q_2_s_i_E1[i,l] = m_2_i_E1[i] * cp_l (Tsat2_E1[l]) * Tsat2_E1[l]/3600
            
#                 #q_2_i[i,l] = q_2_s_i[i,l]
            
#             # Calor transferido a la pared (estacionario)
            
#                 #q_1_s_w[j] = q_2_e_p[l]
    


####################################################
## PASTEURIZADOR
####################################################


def KPI_pasteur(m_suero,
                T_in_suero,
                T_out_suero,
                Tsat_v,
                T_cond,
                f1):
    """
    Calcula:
      - Q_suero [W] integrando cp_suero en función de T
      - m_vapor (kg/h) requerido
      - ratio kg vapor / kg suero

    Parámetros:
    -----------
    m_suero_h : float
        Caudal másico de suero en kg/h
    T_in_suero : float
        Temperatura inicial del suero (°C)
    T_out_suero : float
        Temperatura final deseada (°C)
    f1_solid : float
        Fracción de sólidos totales (0-1)
    Delta_h_vapor_kJ_kg : float
        Energía cedida por 1 kg de vapor (kJ/kg)
        (calor latente + sensible del condensado hasta su T de salida)

    Retorna:
    --------
    dict con:
        'Q_suero': [W] necesarios para calentar el suero
        'm_vapor_req': kg/h de vapor requerido
        'ratio_vapor_producto': kg vapor / kg suero
    """
    # 1. Energía para calentar el suero (kJ/h) integrando cp_suero
    Q_suero = Q_suero_integral(m_suero, T_in_suero, T_out_suero, f1)

    # 2. Masa de vapor requerida
    #    m_vapor = Q_suero / (Delta_h_vapor)
    m_vapor_req = (Q_suero / Delta_h_subenf(Tsat_v,T_cond))*3600

    # 3. Razón vapor/producto
    ratio_vapor_producto = m_vapor_req / m_suero

    return {
        'Q_suero en [kW]': Q_suero/1000,
        'm_vapor_req en [kg/h]': m_vapor_req,
        'ratio_vapor_producto': ratio_vapor_producto
    }

##########################################################
## SECADO
##########################################################

def secado(m_air,T1,T2,T3,T4):
    """Calcula el flujo de vapor requerido para calentar el aire de T1 a T2.

    Parameters
    ----------
    T1 : float
        Temperatura de entrada del aire [°C]
    T2 : float  
        Temperatura de salida del aire en [°C]
    T3 : float          
        Temperatura de entrada vapor en [°C]
    T4 : float      
        Temperatura de salida del condensado en [°C]
    """
    T1_k = T1 + 273
    T2_k = T2 + 273
    P= 101325
    Q_air = (
    m_air * cp_air((T1_k + T2_k) / 2, P) * (T2_k - T1_k) / 3600)  # calor necesario al vapor en [W]

    h_fg = Delta_h_subenf(T3, T4)

    m_req_sec = Q_air / h_fg * 3600

    return m_req_sec


# def pasteurizador_ode(t,T,params):
#     """Genera la ecuación diferencial ordinaria para la temperatura de pasteurización.

#     Parameters
#     ----------
#     t : float
#         Variable independiente (tiempo)
#     T : float
#         Vector de estado
#     params : dict
#         diccionario 
#     """
#     Tpast = T[0]
    
#     # Llamamos a los parámetros desde el diccionario

#     alpha       = params["alpha"]        # fracción de condensación (0 a 1)
#     m_suero     = params["m_suero"]      # flujo de suero (kg/s)
#     cp_suero    = params["cp_suero"]     # calor específico del suero
#     m_hold      = params["m_hold"]       # masa (o volumen*densidad) de suero en el equipo
#     U           = params["U"]            # coef. global de transferencia
#     A           = params["A"]            # área de intercambio
#     T_in        = params["T_in"]         # temp. de entrada del suero
#     T_vapor     = params["T_vapor"]      # temp. del vapor

#     # Cálculo de calor por vapor condensado a través de la fracción alpha
#     Q_vap = alpha*U*A*(T_vapor-Tpast)

#     # Flujo neto de energía por entrada y salida de suero
#     Q_suero = m_suero* cp_suero * (T_in - Tpast)

#     # Balance de energía
#     dTpast_dt = (Q_suero+Q_vap)/ (m_hold * cp_suero)

#     return [dTpast_dt]

# def simulate_pasteurization():
#     # 1) Sorteamos valores "incógnita" según las distribuciones pedidas:
#     T_in_rand     = np.random.uniform(55, 65)   # entre 55 y 65 °C
#     T_vapor_rand  = np.random.uniform(77, 83)   # entre 77 y 83 °C
#     alpha_rand    = np.random.uniform(0, 1)     # fracción (0 a 1)
#     cp_suero_rand = cp_suero(T_in_rand,0.115)
#     m_suero_rand = np.random.uniform(19000,21000)/3600
#     # 2) Definimos los parámetros fijos (¡ajusta según tu sistema real!)
#     params = {
#         "alpha"    : 1, 
#         "m_suero"  : m_suero_rand,      # kg/s, ejemplo
#         "cp_suero" : cp_suero_rand,      # J/(kg·°C) o (J/kgK), revisa la unidad
#         "m_hold"   : 103,    # kg de suero "en el interior" (ejemplo)
#         "U"        : 8000,    # W/(m^2·°C) (?)
#         "A"        : 4.15,      # m^2 de intercambio
#         "T_in"     : T_in_rand,
#         "T_vapor"  : T_vapor_rand
#     }
    
#     # 3) Condición inicial (ejemplo: al inicio todo el pasteurizador está a T_in)
#     T0 = [params["T_in"]]
    
#     # 4) Rango de integración en tiempo (segundos, por ejemplo)
#     t_span = (0, 2000)  # 0 a 2000 s
    
#     # 5) Resolvemos la EDO con solve_ivp
#     sol = solve_ivp(
#         fun=pasteurizador_ode,
#         t_span=t_span,
#         y0=T0,
#         args=(params,),
#         dense_output=True
#     )
    
#     # 6) Evaluamos la solución para graficar
#     t_eval = np.linspace(t_span[0], t_span[1], 300)
#     T_past = sol.sol(t_eval)[0]  # La primera (y única) ecuación
    
#     # 7) Mostramos en pantalla resultados básicos
#     print("=== Parámetros sorteados ===")
#     print(f"T_in       = {params['T_in']:.2f} °C")
#     print(f"T_vapor    = {params['T_vapor']:.2f} °C")
#     print(f"alpha      = {params['alpha']:.2f}")
#     print(f"m_suero    = {params['m_suero']:.2f} kg/s")
#     print(f"cp_suero    = {params['cp_suero']:.2f} j/kgK")
#     print("Solución de la EDO:", "OK" if sol.success else "FALLA")
    
#     # 8) Graficamos la evolución de T_past en el tiempo
#     plt.figure()
#     plt.plot(t_eval, T_past, label='T_past (°C)')
#     plt.axhline(72, color='r', linestyle='--', label='72°C objetivo')
#     plt.xlabel("Tiempo (s)")
#     plt.ylabel("Temperatura (°C)")
#     plt.title("Evolución de la temperatura en la pasteurización")
#     plt.legend()
#     plt.grid(True)
#     plt.show()
#     return sol
