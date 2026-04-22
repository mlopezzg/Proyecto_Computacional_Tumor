# Simulación de Crecimiento de Gliomas mediante Métodos Numéricos

Este repositorio contiene el desarrollo de un modelo computacional para simular la evolución espacial y temporal de tumores cerebrales (gliomas). El proyecto se centra en la implementación de ecuaciones de reacción-difusión para modelar la proliferación y la invasión celular.

## Estructura del Proyecto
* **notebooks/**: Contiene el análisis detallado y las simulaciones en 1D y 2D.
* **codigos/**: Funciones auxiliares y scripts de procesamiento.
* **informes/**: Reportes técnicos con los avances del proyecto.
* **graficas/**: Visualización de los perfiles de densidad tumoral.

## Aspectos Técnicos
El modelo resuelve la ecuación fundamental de crecimiento tumoral:
$$\frac{\partial u}{\partial u} = D \nabla^2 u + \rho u (1 - u)$$

Donde:
* $u$: Densidad de células tumorales.
* $D$: Coeficiente de difusión (invasividad).
* $\rho$: Tasa de proliferación (crecimiento logístico).

Para la resolución numérica se han implementado métodos de diferencias finitas, analizando la estabilidad y convergencia en mallas unidimensionales y bidimensionales.

## Requisitos
* Python 3.x
* NumPy, Matplotlib, SciPy
