## Método Iterativo Gauss-Seidel

Para resolver **Ax = b**, Gauss-Seidel actualiza cada componente secuencialmente:

$$x_i^{\text{new}} = \frac{1}{a_{ii}} \left( b_i - \sum_{j<i} a_{ij} x_j^{\text{new}} - \sum_{j>i} a_{ij} x_j^{\text{old}} \right)$$

### Pseudocódigo

```
x = 0
repeat
    max_diff = 0
    for i = 1..N:
        sigma = sum_{j != i} a[i,j] * x[j]
        x_new = (b[i] - sigma)/a[i,i]
        max_diff = max(max_diff, |x_new - x[i]|)
        x[i] = x_new
until max_diff < tol
```


### Convergencia

La principal diferencia con el método de Jacobi es que Gauss-Seidel utiliza valores recién calculados tan pronto como están disponibles; esta diferencia hace que Gauss-Seidel converja más rápido que Jacobi en la mayoría de los casos.

### Condiciones de convergencia
1. Matriz Diagonalmente Dominante:

$$|a_{ii}| > \sum_{j \neq i} |a_{ij}| \quad \text{para todo } i$$

2. Matriz Simétrica Definida Positiva (SPD):

    - $A = A^T$ (simétrica)
    - Todos los valores propios > 0

3. Matriz con Elementos Positivos en la Diagonal:
Si $a_{ii} > 0$ y la matriz es irreducible, converge.

Propiamente, para el sistema $(I - kA)u^{n+1} = u^n$:

- La matriz $I - kA$ es diagonalmente dominante cuando $k$ es pequeño
- Es simétrica definida positiva porque $A$ (Laplaciano discreto) es SPD
- Esto garantiza la convergencia de Gauss-Seidel

### Ejemplo

Supongamos el sistema:

$$\begin{bmatrix} 4 & 1 \\ 1 & 3 \end{bmatrix} \begin{bmatrix} x_1 \\ x_2 \end{bmatrix} = \begin{bmatrix} 1 \\ 2 \end{bmatrix}$$

Gauss-Seidel actualiza primero $x_1$ usando $x_2$ antiguo, luego $x_2$ usando el nuevo $x_1$. Tras varias iteraciones, converge a la solución exacta.

## Esquema Red-Black

Para mejorar la paralelización y la convergencia, se colorea la malla como un tablero de ajedrez:

- Nodos rojos: $(i + j) \mod 2 = 0$
- Nodos negros: $(i + j) \mod 2 = 1$

Primero se actualizan los nodos rojos usando valores negros antiguos, luego los negros usando valores rojos nuevos. Esto acelera la convergencia respecto al Gauss-Seidel estándar.

## Condiciones de frontera

El código permite elegir condiciones de frontera:

- **Dirichlet homogénea ($u = 0$):** paredes frías, el calor se fuga y la solución tiende a cero.
- **Dirichlet fija ($u = T_{\text{borde}}$):** se fuerzan temperaturas de borde, la solución alcanza equilibrio estacionario.
- **Neumann ($\frac{\partial u}{\partial n} = 0$):** frontera aislada, no hay transferencia de calor hacia afuera, la energía total se conserva.

### Ejemplo físico

Si se coloca un pulso de calor inicial en el centro:

- Con **Dirichlet $u = 0$**, el pulso se disipa rápidamente.
- Con **Dirichlet fija**, el sistema alcanza un estado estacionario con temperaturas de borde constantes.
- Con **Neumann**, el calor permanece dentro del dominio y se redistribuye hasta homogenizarse.


## Flujo completo del programa

1. Se imprimen opciones de frontera y el usuario elige una.
2. Se construye el objeto solver correspondiente.
3. Se genera la condición inicial (ejemplo: distribución gaussiana).
4. Para cada paso de tiempo:
    (a) Se ejecuta el solver iterativo.
    (b) Se guarda un snapshot si corresponde.
    (c) Se analiza máximo y mínimo de la temperatura.

## Relación entre Euler Implícito y Gauss-Seidel

- **Euler implícito** convierte la PDE en un sistema lineal.
- **Gauss-Seidel** (o Red-Black) es el método usado para resolver ese sistema.
- En cada paso temporal se ejecutan múltiples iteraciones de Gauss-Seidel hasta convergencia.

