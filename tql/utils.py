####################################################
## CARGAR DATOS
####################################################

import pandas as pd
import numpy as np

def cargar_datos():
    """Carga datos de flujo másico de vapor en [kg/h] y genera flujo másico de producto en [kg/h].

    Returns
    -------
    (tuple)
        Cada elemento de la tupla es un array y queda definido como a,b = cargar_datos()
    Ejemplo como cargar: a,b = cargar_datos()
    con a =  mv_evap y b = mprd_in
    """
    mv_evap = pd.read_csv('Mvap_i.csv')/3
    mv_evap = mv_evap['Mvap'].to_numpy()
    mprd_in = np.random.uniform(19000,21000,len(mv_evap))
    return mv_evap,mprd_in