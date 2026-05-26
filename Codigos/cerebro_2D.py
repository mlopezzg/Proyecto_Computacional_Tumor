import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np
import os

# 1. Cargar el archivo MINC
directorio_script = os.path.dirname(os.path.abspath(__file__))

archivo_mnc = os.path.join(directorio_script, 't1_icbm_normal_1mm_pn0_rf0.mnc.gz')

img = nib.load(archivo_mnc)

# 2. Extraer la matriz de datos numéricos (array 3D)
data = img.get_fdata()

# 3. Extraer un corte axial (plano XY). En este MINC, el eje axial es el índice 0.
corte_z = data.shape[0] // 2
corte_axial = data[corte_z, :, :]


''' 
# BLOQUE GRÁFICO 1 SILENCIADO
# 4. Visualizar la matriz en escala de grises
plt.figure(figsize=(6, 6))
plt.imshow(corte_axial, cmap='gray')
plt.title(f'Corte Axial del Cerebro (Z = {corte_z})')
plt.colorbar(label='Intensidad de señal (T1)')
plt.axis('off')
plt.show()
'''

# --- NUEVA SECCIÓN: ANÁLISIS DE INTENSIDADES ---

# Aplanamos la matriz 2D a un vector 1D para hacer el histograma.
# Filtramos los valores muy bajos (> 10) para ignorar el fondo negro gigante 
# y enfocarnos solo en el cerebro.
datos_pixeles = corte_axial[corte_axial > 10].flatten()

'''
# BLOQUE GRÁFICO 2 SILENCIADO
# Creamos el histograma
plt.figure(figsize=(8, 5))
plt.hist(datos_pixeles, bins=100, color='gray', edgecolor='black', alpha=0.7)
plt.title('Histograma de Intensidades del Corte Axial (T1)')
plt.xlabel('Valor de Intensidad del Píxel')
plt.ylabel('Frecuencia (Cantidad de píxeles)')
plt.grid(True, alpha=0.5)
plt.show()
'''

# --- NUEVA SECCIÓN: CREACIÓN DE LA MATRIZ DE DIFUSIÓN D(x,y) ---

# 1. Definimos las velocidades de difusión biológicas reales (ejemplo)
# La célula tumoral migra ~5 veces más rápido en sustancia blanca.
D_gris = 0.01
D_blanca = 0.05
D_cero = 0.0

# 2. Creamos una matriz vacía del mismo tamaño que el cerebro
matriz_D = np.zeros_like(corte_axial, dtype=float)

# 3. Aplicamos las reglas de umbrales según el histograma
for i in range(corte_axial.shape[0]):
    for j in range(corte_axial.shape[1]):
        pixel = corte_axial[i, j]
        
        if 150 <= pixel <= 400:
            matriz_D[i, j] = D_gris
        elif pixel > 600:
            matriz_D[i, j] = D_blanca
        else:
            matriz_D[i, j] = D_cero


# BLOQUE GRÁFICO 3 SILENCIADO
# 4. Visualizamos la matriz D resultante (para verificar la segmentación)
plt.figure(figsize=(6, 6))
plt.imshow(matriz_D, cmap='viridis') # 'viridis' resalta los 3 niveles de D
plt.title('Matriz de Difusión Espacial D(x,y)')
plt.colorbar(label='Coeficiente de Difusión D')
plt.axis('off')
plt.show()
