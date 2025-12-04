# Metodología Numérica

## Esquema de Euler Implícito

Para resolver la ecuación de calor numéricamente, utilizamos el **método de diferencias finitas** con el **esquema de Euler Implícito**.

## Discretización del Dominio

### Malla Espacial

Dividimos el dominio \([0, a] \times [0, b]\) en una malla uniforme de \(N \times N\) nodos interiores. Se definen los siguientes pasos espaciales:

\[
\begin{align}
h_x = \frac{a}{N+1} \\
h_y = \frac{b}{N+1}
\end{align}
\]

### Malla Temporal

El intervalo temporal \([0, T]\) se divide en \(M\) pasos donde \(k = \frac{T}{M}\) es el paso temporal.

## Aproximación por Diferencias Finitas

### Derivadas Espaciales (Diferencias Centrales)

Para las segundas derivadas espaciales usamos diferencias centrales de segundo orden:

\[
\frac{\partial^2 u}{\partial x^2}\bigg|_{i,j} \approx \frac{u_{i-1,j} - 2u_{i,j} + u_{i+1,j}}{h_x^2}
\]

\[
\frac{\partial^2 u}{\partial y^2}\bigg|_{i,j} \approx \frac{u_{i,j-1} - 2u_{i,j} + u_{i,j+1}}{h_y^2}
\]

**Error de truncamiento**: \(\mathcal{O}(h^2)\)

### Derivada Temporal (Diferencias hacia adelante)

Se aproxima la derivada temporal haciendo diferencias hacia adelante:

\[
\frac{\partial u}{\partial t}\bigg|_{i,j}^{n} \approx \frac{u_{i,j}^{n+1} - u_{i,j}^n}{k}
\]


## Esquema Numérico Completo

### Ecuación Discretizada

Se sustituyen las aproximaciones anteriores en la ecuacion de calor:

\[
\frac{u_{i,j}^{n+1} - u_{i,j}^n}{k} = c^2 \left[\frac{u_{i-1,j}^{n+1} - 2u_{i,j}^{n+1} + u_{i+1,j}^{n+1}}{h_x^2} + \frac{u_{i,j-1}^{n+1} - 2u_{i,j}^{n+1} + u_{i,j+1}^{n+1}}{h_y^2}\right]
\]

De lo anterior, se puede definir:

\[
\mu_x = \frac{c^2 k}{h_x^2}, \quad \mu_y = \frac{c^2 k}{h_y^2}
\]

Reordenando queda:

\[
(1 + 2\mu_x + 2\mu_y)u_{i,j}^{n+1} - \mu_x(u_{i-1,j}^{n+1} + u_{i+1,j}^{n+1}) - \mu_y(u_{i,j-1}^{n+1} + u_{i,j+1}^{n+1}) = u_{i,j}^n
\]

donde el lado derecho (RHS) contiene solo términos conocidos en el tiempo \(t_n\).

## Sistema Lineal

### Estructura Matricial

El sistema completo obtenido anteriormente se puede escribir como:

\[
A \mathbf{u}^{n+1} = \mathbf{b}^n
\]

donde:

- \(A\): matriz de coeficientes.
- \(\mathbf{u}^{n+1}\): vector de incógnitas en el tiempo \(t_{n+1}\).
- \(\mathbf{b}^n\): vector del lado derecho (conocido).

### Ordenamiento de Incógnitas

Se ordenan las incognitas de la siguiente manera:

\[
\mathbf{u^{n+1}} = [u_{1,1}, u_{2,1}, \ldots, u_{N_x,1}, u_{1,2}, u_{2,2}, \ldots, u_{N_x,2}, \ldots, u_{N_x,N_y}]^T
\]

Este ordenamiento resulta en una matriz \(A\) de dimensión \((N \times N) \times (N \times N)\).

### Estructura de la Matriz \(A\)

De lo obtenido anteriormente, podemos ver que para un punto $(i,j)$ dentro de la malla, su ecuacion tambien involucra a sus vecinos horizontales y verticales, donde podemos ver que sus coeficientes se ven de la siguiente manera:

- \(u^{n+1}_{i,j}: 1 + 2\mu_x + 2\mu_y\)
- \({u}^{n+1}_{i\pm1,j}:-\mu_x\)
- \({u}^{n+1}_{i,j\pm1}:-\mu_y\)

Entonces, $A$ es una matriz tridiagonal por bloques con la siguiente estructura:

#### Bloque diagonal B
Corresponden a los nodos de una misma fila, donde su diagonal principal es \(1 + 2\mu_x + 2\mu_y\) y subdiagonal y superdiagonal es \(-\mu_x\). Es decir: \(B= tridiagonal(-\mu_x, 1 + 2\mu_x + 2\mu_y, -\mu_x)\)

#### Bloques fuera de la diagonal
Son los que conectan filas adyacentes, o sea \(j\) y \(j\pm1\). Entonces cada nodo \((i,j)\) se conecta solo con nodos de la misma columna pero diferente fila. Por lo tanto, estos bloques estan compuestos por matrices diagonales con \(-\mu_y\) en la diagonal, es decir, \(C=-\mu_y I\)

Entonces, la matriz A tiene la siguiente forma:

\[
A= \begin{bmatrix}    
B & C & 0 & \ldots & 0 \\
C & B & C & \ldots & 0 \\
0 & C & B & \ddots & \vdots \\
\vdots & \ddots & \ddots & \ddots & C \\    
0 & \ldots & 0 & C & B    
\end{bmatrix}    
\]

Esta es una **matriz dispersa** que se puede almacenar y resolver eficientemente.
