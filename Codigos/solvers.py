import numpy as np

def ftcs_1d(D, rho, nx, dx, dt, nt, u_inicial):
    """
    Resuelve la ecuación de Fisher-KPP en 1D mediante Diferencias Finitas Explícitas (FTCS).
    Retorna la matriz final de densidades y las listas de seguimiento cinemático.
    """
    u = u_inicial.copy()
    tiempos = []
    posiciones_frente = []

    for n in range(nt):
        un = u.copy()
        
        # Operador Laplaciano 1D vectorizado y término de reacción logística
        difusion = D * (un[2:] - 2*un[1:-1] + un[:-2]) / dx**2
        reaccion = rho * un[1:-1] * (1 - un[1:-1])
        
        u[1:-1] = un[1:-1] + dt * (difusion + reaccion)
        
        # Seguimiento de la cinemática del frente de onda
        indices_tumor = np.where(u > 0.5)[0]
        if len(indices_tumor) > 0:
            frente_x = indices_tumor[-1] * dx
            tiempo_actual = n * dt
            # Registro de datos tras superar el transiente inicial
            if tiempo_actual > 15.0: 
                tiempos.append(tiempo_actual)
                posiciones_frente.append(frente_x)

    return u, tiempos, posiciones_frente


def ftcs_2d_homogeneo(D, rho, nx, ny, dx, dy, dt, nt, u_inicial):
    """
    Resuelve la ecuación de Fisher-KPP en 2D para un medio homogéneo
    usando Diferencias Finitas Explícitas (FTCS).
    Retorna la matriz bidimensional con la densidad celular en el tiempo final.
    """
    u = u_inicial.copy()
    
    for n in range(nt):
        un = u.copy()
        
        # Operador Laplaciano 2D vectorizado
        difusion_x = (un[2:, 1:-1] - 2*un[1:-1, 1:-1] + un[:-2, 1:-1]) / dx**2
        difusion_y = (un[1:-1, 2:] - 2*un[1:-1, 1:-1] + un[1:-1, :-2]) / dy**2
        
        laplaciano = difusion_x + difusion_y
        
        # Término de reacción logística
        reaccion = rho * un[1:-1, 1:-1] * (1 - un[1:-1, 1:-1])
        
        # Actualización explícita
        u[1:-1, 1:-1] = un[1:-1, 1:-1] + dt * (D * laplaciano + reaccion)

    return u