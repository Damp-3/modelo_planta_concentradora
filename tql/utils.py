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



def resample_xlsx_inst(file_path, resample_interval='5min', output_file='resampled_output.csv'):
    """
    Función para leer un archivo .xlsx con datos de tiempo y valores de medición,
    resamplear a un intervalo fijo y realizar interpolación lineal.
    
    Parámetros:
    - file_path (str): Ruta del archivo .xlsx de entrada.
    - resample_interval (str): Intervalo de resampleo (ejemplo: '5min').
    - output_file (str): Nombre del archivo de salida con los datos resampleados en formato CSV.
    
    Retorna:
    - DataFrame con los datos resampleados.
    """
    # Cargar datos
    df = pd.read_excel(file_path, dtype={'Fecha': str})
    
    # Convertir la columna de fecha a formato datetime
    df['Fecha'] = pd.to_datetime(df['Fecha'])
    
    # Renombrar columna de valores para consistencia
    df.columns = ['Fecha', 'Valor']
    
    # Resamplear la serie temporal con el intervalo especificado
    df_resampled = df.set_index('Fecha').resample(resample_interval).mean()
    
    # Interpolación lineal para rellenar valores faltantes
    df_resampled = df_resampled.interpolate(method='linear', limit_direction='both')
    
    # Eliminar la columna de fecha, dejando solo el índice y valores
    df_resampled.reset_index(drop=True, inplace=True)
    
    # Guardar el archivo de salida en formato CSV
    df_resampled.to_csv(output_file, index=False)
    

    return df_resampled

def resample_xlsx_hist(file_path, resample_interval='5min', output_file='resampled_output.csv'):
    """
    Función para leer un archivo .xlsx con datos de tiempo y valores de medición,
    resamplear a un intervalo fijo y realizar interpolación lineal.
    
    Parámetros:
    - file_path (str): Ruta del archivo .xlsx de entrada.
    - resample_interval (str): Intervalo de resampleo (ejemplo: '5min').
    - output_file (str): Nombre del archivo de salida con los datos resampleados en formato CSV.
    
    Retorna:
    - DataFrame con los datos resampleados.
    """
    # Cargar datos
    df = pd.read_excel(file_path, dtype={'Fecha': str})
    
    # Convertir la columna de fecha a formato datetime
    df['Fecha'] = pd.to_datetime(df['Fecha'])
    
    # Renombrar columna de valores para consistencia
    df.columns = ['Fecha', 'Valor']
    
    # Resamplear la serie temporal con el intervalo especificado
    df_resampled = df.set_index('Fecha').resample(resample_interval).mean()
    
    # Interpolación lineal para rellenar valores faltantes
    df_resampled = df_resampled.interpolate(method='linear', limit_direction='both')
    
    # Restamos el valor anterior, para convertirlo en un valor instantáneo
    df_resampled = df_resampled.diff().fillna(0)

    # Eliminar la columna de fecha, dejando solo el índice y valores
    df_resampled.reset_index(drop=True, inplace=True)
    
    # Guardar el archivo de salida en formato CSV
    df_resampled.to_csv(output_file, index=False)
    
    return df_resampled

