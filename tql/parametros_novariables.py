##########################################
## PARAMETROS NO VARIABLES
##########################################

import math

def cargar_parametros_novariables():
    """
    Crea un diccionario con los parámetros no variables.

    Returns
    -------
    (dict)
        Diccionario con todos los parámetros no variables.

    Examples
    -------

    parametros = cargar_parametros_novariables()
    print(parametros)

    """
    # Constantes fisicas
    g = 9.81

    # Constantes geométricas
    do = 0.0508  # Diámetro externo [m]
    e = 0.005    # Espesor [m]
    di = do - 2 * e  # Diámetro interno [m]
    
    # Configuración de tubos
    tubos_por_efecto = {
        "1er efecto": 180,
        "2do efecto": 180,
        "3er efecto": 120,
        "4to efecto": 120
    }
    
    # Longitud de los tubos
    L = 13  # Longitud [m]
    
    # Área intercambio pasteurizador
    A_pasteur = math.pi*L*do*2

    # Masa en el pasteurizador
    m_hold = math.pi*(di**2)*L*1060

    # Coeficiente global de transferencia de calor pasteurizador
    U = 8000 

    # Área transversal y longitudinal de tubos
    at_un_tubo = math.pi * (do ** 2) / 4  # Área transversal de un tubo [m²]
    al_un_tubo = math.pi * do * L        # Área longitudinal de un tubo [m²]
    
    # Configuración de evaporadores
    diametros_evaporadores = {
        "1er efecto": 1.00,
        "2do efecto": 1.00,
        "3er efecto": 0.95,
        "4to efecto": 0.95
    }
    # Temperaturas saturacion vapor
    tsat1 = {
        "1er efecto":353,
        "2do efecto":343,
        "3er efecto":333,
        "4to efecto":323
    }
    # Cálculos
    parametros = {
        "a_p": A_pasteur,
        "U_past":U,
        "m_hold":m_hold,
        "tsat1":tsat1,
        "L": L,
        "g": g,
        "di": di,
        "at_un_tubo": at_un_tubo,
        "al_un_tubo": al_un_tubo,
        "tubos_por_efecto": tubos_por_efecto,
        "areas_transversales_totales": {
            efecto: at_un_tubo * n_tubos
            for efecto, n_tubos in tubos_por_efecto.items()
        },
        "areas_longitudinales_totales": {
            efecto: al_un_tubo * n_tubos
            for efecto, n_tubos in tubos_por_efecto.items()
        },
        "areas_evaporadores": {
            efecto: {
                "area_transversal": math.pi * (d ** 2) / 4,
                "area_longitudinal": math.pi * d * L
            }
            for efecto, d in diametros_evaporadores.items()
        }
    }
    
    return parametros
