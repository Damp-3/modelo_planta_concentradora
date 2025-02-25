#########################################
##### VALIDATION DATA
#########################################

import numpy as np
import pandas as pd
from tql.utils import *
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from scipy.signal import correlate
from statsmodels.tsa.stattools import grangercausalitytests

# Análisis de correlación lineal
# Coeficiente de Pearson

def pearson_correlation_matrix(*dataframes):
    """
    Calcula la matriz de correlación de Pearson para múltiples DataFrames y genera un heatmap.
    
    Parámetros:
    - *dataframes: Múltiples DataFrames con una única columna de datos.
    
    Retorna:
    - DataFrame con la matriz de correlación de Pearson.
    - Muestra un heatmap con las correlaciones.
    """
    # Convertir los DataFrames en un solo DataFrame con nombres de columnas únicos
    df_combined = pd.concat(dataframes, axis=1)
    
    # Asignar nombres de variables si no tienen nombres definidos
    df_combined.columns = ['temp agua','q gas','q agua','m vap','t vap','p vap']
    
    # Calcular la matriz de correlación de Pearson
    correlation_matrix = df_combined.corr()
    
    # Graficar el heatmap
    plt.figure(figsize=(8,6))
    sns.heatmap(correlation_matrix, annot=True, cmap="coolwarm", fmt=".2f", linewidths=0.5)
    plt.title("Matriz de Correlación de Pearson")
    plt.show()
    
    return correlation_matrix


# Correlación cruzada

def cross_correlation_analysis(array1, array2, max_lag=50):
    """
    Calcula la correlación cruzada entre dos arrays y analiza el desfase temporal.
    
    Parámetros:
    - array1: Primer array de datos.
    - array2: Segundo array de datos.
    - max_lag: Máximo desfase en ambas direcciones.
    
    Retorna:
    - Lags: Desfase temporal en índices.
    - Cross-correlation values: Valores de correlación cruzada.
    - Muestra un gráfico de correlación cruzada.
    """
    # Asegurar que los arrays sean 1D
    array1 = np.asarray(array1).flatten()
    array2 = np.asarray(array2).flatten()
    
    # Verificar que ambos arrays tengan la misma longitud
    if len(array1) != len(array2):
        raise ValueError("Los arrays deben tener la misma longitud.")
    
    # Normalizar los datos (opcional, pero recomendado para mejorar la interpretación)
    array1 = (array1 - np.mean(array1)) / np.std(array1)
    array2 = (array2 - np.mean(array2)) / np.std(array2)
    
    # Calcular la correlación cruzada
    cross_corr = correlate(array1, array2, mode='full')
    lags = np.arange(-len(array1) + 1, len(array1))
    
    # Limitar a max_lag
    mid = len(lags) // 2
    lags = lags[mid - max_lag:mid + max_lag + 1]
    cross_corr = cross_corr[mid - max_lag:mid + max_lag + 1]
    
    # Graficar correlación cruzada
    plt.figure(figsize=(8, 4))
    plt.plot(lags, cross_corr, marker='o')
    plt.axvline(0, color='k', linestyle='--', label='Lag 0')
    plt.xlabel("Desfase (lags)")
    plt.ylabel("Correlación cruzada")
    plt.title("Análisis de Correlación Cruzada")
    plt.legend()
    plt.grid()
    plt.show()
    
    return lags, cross_corr

# Granger Causality Analisis

def granger_causality_analysis(array1, array2, max_lag=5, significance_level=0.05):
    """
    Realiza el test de causalidad de Granger entre dos series temporales.
    
    Parámetros:
    - array1: Primer array de datos (posible causa).
    - array2: Segundo array de datos (posible efecto).
    - max_lag: Máximo número de retardos a considerar.
    - significance_level: Nivel de significancia para determinar causalidad.
    
    Retorna:
    - Resultados del test de Granger.
    """
    # Asegurar que los arrays sean 1D
    array1 = np.asarray(array1).flatten()
    array2 = np.asarray(array2).flatten()
    
    # Verificar que ambos arrays tengan la misma longitud
    if len(array1) != len(array2):
        raise ValueError("Los arrays deben tener la misma longitud.")
    
    # Convertir a DataFrame
    df = pd.DataFrame({'Variable1': array1, 'Variable2': array2})
    
    # Realizar el test de causalidad de Granger
    print(f"Ejecutando test de Granger con max_lag={max_lag} y nivel de significancia={significance_level}")
    test_results = grangercausalitytests(df, max_lag, verbose=True)
    
    return test_results


