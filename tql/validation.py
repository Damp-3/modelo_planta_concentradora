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
from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_absolute_error, mean_squared_error

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
    cross_corr = correlate(array1, array2, mode='full')/len(array1)
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


def linear_regression_feature_importance(X, y):
    """
    Realiza una regresión lineal múltiple y muestra la importancia de las variables.
    """
    model = LinearRegression()
    model.fit(X, y)
    importance = model.coef_
    plt.figure(figsize=(8,6))
    plt.bar(range(len(importance)), importance)
    plt.xlabel("Variables")
    plt.ylabel("Importancia de coeficientes")
    plt.title("Importancia de Variables en Regresión Lineal")
    plt.show()
    return importance

def random_forest_feature_importance(X, y):
    """
    Entrena un modelo Random Forest y muestra la importancia de las características.
    """
    model = RandomForestRegressor(n_estimators=100, random_state=42)
    model.fit(X, y)
    importance = model.feature_importances_
    plt.figure(figsize=(8,6))
    plt.bar(range(len(importance)), importance)
    plt.xlabel("Variables")
    plt.ylabel("Importancia")
    plt.title("Importancia de Variables en Random Forest")
    plt.show()
    return importance

import numpy as np
import matplotlib.pyplot as plt

def plot_multiple_series_normalized(arrays, labels, title="Comparación de Series de Datos Normalizadas"):
    """
    Plotea múltiples series de datos normalizadas en un solo gráfico.
    
    Parámetros:
    - arrays: Lista de arrays con los datos a graficar.
    - labels: Lista de etiquetas para cada serie.
    - title: Título del gráfico.
    """
    if len(arrays) != len(labels):
        raise ValueError("El número de etiquetas debe coincidir con el número de series de datos.")
    
    plt.figure(figsize=(10, 5))
    
    for i, data in enumerate(arrays):
        normalized_data = data / np.max(data)  # Normalización por el valor máximo
        plt.plot(normalized_data, label=labels[i])
    
    plt.xlabel("Tiempo")
    plt.ylabel("Valor Normalizado")
    plt.title(title)
    plt.legend()
    plt.grid()
    plt.show()

def compare_resampled_interpolation(csv_file, time_column, value_column, resample_interval='5min'):
    """
    Compara los datos originales con los datos resampleados e interpolados.

    Parámetros:
    - csv_file: Ruta del archivo CSV.
    - time_column: Nombre de la columna de tiempo.
    - value_column: Nombre de la columna de valores a analizar.
    - resample_interval: Intervalo de resampleo (por defecto, '5min').

    Retorna:
    - Gráfico con datos originales vs resampleados/interpolados.
    """
    # Cargar datos
    df = pd.read_excel(csv_file)

    df[time_column] = pd.to_datetime(df[time_column])
    df = df.sort_values(by=time_column).set_index(time_column)
    
    df_resampled = df[[value_column]].resample(resample_interval).mean()
    df_resampled[value_column] = df_resampled[value_column].interpolate(method='linear')
    
    # Ajustar los arrays al mismo tamaño sin alinearlos
    min_length = min(len(df[value_column].dropna()), len(df_resampled[value_column].dropna()))
    original_values = df[value_column].dropna().values[:min_length]
    resampled_values = df_resampled[value_column].dropna().values[:min_length]
    
    # Calcular métricas de comparación
    mae = mean_absolute_error(original_values, resampled_values)
    rmse = np.sqrt(mean_squared_error(original_values, resampled_values))
    
    print(f"Error absoluto medio (MAE): {mae:.4f}")
    print(f"Raíz del error cuadrático medio (RMSE): {rmse:.4f}")
    
    plt.figure(figsize=(10, 5))
    plt.plot(df.index[:min_length], original_values, label="Datos Originales", linestyle='-', marker='o', alpha=0.6)
    plt.plot(df_resampled.index[:min_length], resampled_values, label="Datos Resampleados + Interpolados", linestyle='--', marker='x', alpha=0.8)
    
    plt.xlabel("Tiempo")
    plt.ylabel(value_column)
    plt.title(f"Comparación de {value_column}: Original vs Resampleado + Interpolado")
    plt.legend()
    plt.grid()
    plt.show()
    
    return mae, rmse