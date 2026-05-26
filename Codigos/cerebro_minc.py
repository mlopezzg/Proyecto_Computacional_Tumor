import nibabel as nib
import matplotlib.pyplot as plt
import numpy as np

# 1. Cargar el archivo MINC
archivo_mnc = 'Codigos/t1_icbm_normal_1mm_pn0_rf0.mnc.gz'
img = nib.load(archivo_mnc)

# 2. Extraer la matriz de datos numéricos (array 3D)
data = img.get_fdata()

# 3. Extraer un corte axial (plano XY). En este MINC, el eje axial es el índice 0.
corte_z = data.shape[0] // 2
corte_axial = data[corte_z, :, :]



# 5. Visualizar la matriz en escala de grises
plt.figure(figsize=(6, 6))
plt.imshow(corte_axial, cmap='gray')
plt.title(f'Corte Axial del Cerebro (Z = {corte_z})')
plt.colorbar(label='Intensidad de señal (T1)')
plt.axis('off')
plt.show()