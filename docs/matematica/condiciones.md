## Tipos de Condiciones

### 1. Condiciones de Dirichlet

#### Implementación Numérica

Para una frontera izquierda con temperatura \(g\):

\[
u_{0,j}^{n+1} = g(y_j, t_{n+1})
\]

Este es un valor conocido, por lo que se mueve al lado derecho del sistema:

```python
# Punto interior más cercano a la frontera izquierda
idx = j * N + 0  # i = 0

# Término que se añade al RHS
b[idx] += mu_x * g_left(y_j, t_n)
```

La logica es la misma los otros bordes.

### 2. Condiciones de Neumann

## Implementación numérica con Puntos Fantasma

Para implementar condiciones de Neumann, introducimos **puntos fantasma** fuera del dominio.

#### Frontera Izquierda (\(x = 0\))

La condición de Neumann es:

\[
\frac{\partial u}{\partial x}\bigg|_{x=0} = g(y, t)
\]

Se aproxima esta derivada por diferencias centrales:

\[
\frac{u_{1,j} - u_{-1,j}}{2h_x} = g(y_j, t)
\]

donde \(u_{-1,j}\) es el punto fantasma.

Se despeja este punto fantasma::

\[
u_{-1,j} = u_{1,j} - 2h_x \cdot g(y_j, t)
\]

Se sustituye en la ecuacion de calor discretizada:

\[
\begin{align}
(1 + 2\mu_x + 2\mu_y)u_{0,j}^{n+1} - \mu_x(2u_{1,j}^{n+1} - 2h_x g) - \mu_y(u_{0,j-1}^{n+1} + u_{0,j+1}^{n+1}) = u_{0,j}^n
\end{align}
\]

\[
(1 + 2\mu_x + 2\mu_y)u_{0,j}^{n+1} - 2\mu_xu_{1,j}^{n+1}  - \mu_y(u_{0,j-1}^{n+1} + u_{0,j+1}^{n+1}) = u_{0,j}^n - 2h_x g \mu_x
\]

### Modificaciones al Sistema Lineal

Como se puede ver, los resultados anteriores modifican el planteamiento del sistema lineal que obtuvimos en la derivacion anterior:

#### En la Matriz \(A\)

El coeficiente del **vecino** en la dirección de la frontera se **duplica**:

**Antes (Dirichlet o punto interior)**:
```
A[idx, idx]     =  1 + 2μₓ + 2μᵧ
A[idx, idx+1]   = -μₓ           (vecino derecho)
```

**Ahora (Neumann izquierda)**:
```
A[idx, idx]     =  1 + 2μₓ + 2μᵧ
A[idx, idx+1]   = -2μₓ          (coeficiente duplicado)
```

#### En el Vector \(b\) (RHS)

Se añade el término del flujo:

```python
b[idx] -= 2 * mu_x * h_x * g(y_j, t_n)
```

**Signo**: Depende de la orientación de la normal exterior:
- Frontera **izquierda**: \(-2\mu_x h_x g\)
- Frontera **derecha**: \(+2\mu_x h_x g\)
- Frontera **inferior**: \(-2\mu_y h_y g\)
- Frontera **superior**: \(+2\mu_y h_y g\)

## Tabla de Modificaciones

Se presenta la siguiente tabla a modo de resumen:

| Frontera | Normal | Coeficiente duplicado | Término RHS |
|----------|--------|----------------------|-------------|
| Izquierda (\(x=0\)) | \(-\hat{x}\) | \(A[idx, idx+1] = -2\mu_x\) | \(-2\mu_x h_x g\) |
| Derecha (\(x=a\)) | \(+\hat{x}\) | \(A[idx, idx-1] = -2\mu_x\) | \(+2\mu_x h_x g\) |
| Inferior (\(y=0\)) | \(-\hat{y}\) | \(A[idx, idx+N] = -2\mu_y\) | \(-2\mu_y h_y g\) |
| Superior (\(y=b\)) | \(+\hat{y}\) | \(A[idx, idx-N] = -2\mu_y\) | \(+2\mu_y h_y g\) |
