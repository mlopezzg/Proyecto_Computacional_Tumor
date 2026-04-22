import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.widgets import Slider
from matplotlib.gridspec import GridSpec

# PARÁMETROS FÍSICOS Y DE SIMULACIÓN
N = 125                           # Número de partículas
L = 10e-6                         # Dimensiones iniciales de la caja (m)
T0 = 300.0                        # Temperatura inicial (K)
dt = 1e-12                        # Paso de integración (1 ps)
m = 3.32e-27                      # Masa del H2 (kg)
kB = 1.380649e-23                 # Constante de Boltzmann (J/K)

# INICIALIZACIÓN AlEATORIA DE POSICIONES Y VELOCIDADES
pos = np.random.rand(N, 3) * L
vel = np.random.randn(N, 3)

v_rms_deseada = np.sqrt(3 * kB * T0 / m)
v_rms_actual = np.sqrt(np.mean(np.sum(vel**2, axis=1)))
vel = vel * (v_rms_deseada / v_rms_actual)

historial_T = np.full(1000, T0)
historial_P = np.full(1000, (N * kB * T0) / (L**3))
pasos = np.arange(1000)

# INTERFAZ
fig = plt.figure(figsize=(15, 8)) 

plt.subplots_adjust(left=0.05, bottom=0.25, right=0.95, top=0.95, wspace=0.3, hspace=0.4)


gs = GridSpec(2, 3, figure=fig)

# La caja 3D ocupa las dos filas de la primera columna
ax_3d = fig.add_subplot(gs[:, 0], projection='3d')

ax_T = fig.add_subplot(gs[0, 1])
ax_P = fig.add_subplot(gs[0, 2])
ax_E = fig.add_subplot(gs[1, 1])
ax_Hist = fig.add_subplot(gs[1, 2])

# Inicializar gráficos 2D
line_T, = ax_T.plot(pasos, historial_T, 'r-')
line_P, = ax_P.plot(pasos, historial_P, 'b-')
line_E, = ax_E.plot(pasos, historial_T * (1.5 * N * kB), 'g-')

ax_T.set_title('Temperatura T(t)')
ax_P.set_title('Presión P(t)')
ax_E.set_title('Energía Cinética E_k(t)')

# Inicializar gráfico 3D (La Caja)
# alpha=0.6 hace las partículas un poco transparentes para ver las de atrás
scatter_3d = ax_3d.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c='teal', s=30, alpha=0.6)
ax_3d.set_title('Simulación de H2 en Tiempo Real')
ax_3d.set_xlabel('X (m)')
ax_3d.set_ylabel('Y (m)')
ax_3d.set_zlabel('Z (m)')

# SLIDERS
axcolor = 'lightgoldenrodyellow'
ax_slider_L = plt.axes([0.15, 0.10, 0.65, 0.03], facecolor=axcolor)
ax_slider_T = plt.axes([0.15, 0.05, 0.65, 0.03], facecolor=axcolor)

slider_L = Slider(ax_slider_L, 'Lado Caja (μm)', 5.0, 20.0, valinit=L*1e6)
slider_T = Slider(ax_slider_T, 'Termostato (K)', 10.0, 1000.0, valinit=T0)

def update_L(val):
    global L
    L = slider_L.val * 1e-6
slider_L.on_changed(update_L)

def update_T(val):
    global vel
    T_objetivo = slider_T.val
    v_rms_deseada = np.sqrt(3 * kB * T_objetivo / m)
    v_rms_actual = np.sqrt(np.mean(np.sum(vel**2, axis=1)))
    vel = vel * (v_rms_deseada / v_rms_actual)
slider_T.on_changed(update_T)


# Simulación y Animación
def animar(frame):
    global pos, vel, historial_T, historial_P
    
    # ACELERADOR VISUAL: Procesamos 100 pasos físicos de 1 ps por cada cuadro de la animación
    pasos_por_frame = 100
    
    for _ in range(pasos_por_frame):
        # Física de colisiones
        pos += vel * dt
        
        for i in range(3):
            rebote_max = pos[:, i] > L
            rebote_min = pos[:, i] < 0
            
            vel[rebote_max, i] *= -1
            pos[rebote_max, i] = 2*L - pos[rebote_max, i]
            
            vel[rebote_min, i] *= -1
            pos[rebote_min, i] *= -1

   
    v_cuadrado = np.sum(vel**2, axis=1)
    E_k_total = 0.5 * m * np.sum(v_cuadrado)
    T_actual = (2.0 / 3.0) * E_k_total / (N * kB)
    V_actual = L**3
    P_actual = (N * kB * T_actual) / V_actual
    
    historial_T[:-1] = historial_T[1:]
    historial_T[-1] = T_actual
    historial_P[:-1] = historial_P[1:]
    historial_P[-1] = P_actual
    
    # Actualización de gráficos 2D
    line_T.set_ydata(historial_T)
    ax_T.set_ylim(min(historial_T)*0.9, max(historial_T)*1.1)
    
    line_P.set_ydata(historial_P)
    ax_P.set_ylim(min(historial_P)*0.9, max(historial_P)*1.1)
    
    line_E.set_ydata(historial_T * (1.5 * N * kB))
    ax_E.set_ylim(min(historial_T * (1.5 * N * kB))*0.9, max(historial_T * (1.5 * N * kB))*1.1)
    
    ax_Hist.clear()
    ax_Hist.set_title('Distribución de Velocidades')
    magnitud_v = np.sqrt(v_cuadrado)
    ax_Hist.hist(magnitud_v, bins=15, density=True, color='skyblue', edgecolor='black')
    
    v_teo = np.linspace(0, max(magnitud_v)*1.5, 100)
    f_v = (4 * np.pi * v_teo**2) * ((m / (2 * np.pi * kB * T_actual))**1.5) * np.exp((-m * v_teo**2) / (2 * kB * T_actual))
    ax_Hist.plot(v_teo, f_v, 'r-', lw=2)

    # Actualizar la Caja 3D
    scatter_3d._offsets3d = (pos[:, 0], pos[:, 1], pos[:, 2])
    
    ax_3d.set_xlim(0, L)
    ax_3d.set_ylim(0, L)
    ax_3d.set_zlim(0, L)

    return line_T, line_P, line_E, ax_Hist, scatter_3d

# Intervalo un poco más rápido (10 ms) para que se vea fluido
ani = FuncAnimation(fig, animar, frames=200, interval=10, blit=False)

plt.show()