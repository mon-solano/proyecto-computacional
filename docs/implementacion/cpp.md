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

- # Paralelización Red–Black con Bloques (Tiling)

En la resolución numérica de la ecuación de calor en dos dimensiones mediante el método de Euler implícito, cada paso de tiempo requiere resolver un sistema lineal grande. Este sistema se resuelve con Gauss–Seidel, pero su versión clásica es totalmente secuencial, ya que cada punto depende de las actualizaciones recientes de sus vecinos. Esto hace que el método sea lento cuando la malla es grande. Para aprovechar los procesadores multinúcleo y reducir los tiempos de cómputo, se combinan dos técnicas muy utilizadas en paralelización: el esquema Red–Black Gauss–Seidel y la división de la malla en bloques o tiles.

## Método Red–Black

La idea del método Red–Black es dividir la malla como si fuera un tablero de ajedrez: los nodos donde (i + j) es par se marcan como rojos, y los nodos donde (i + j) es impar como negros. Esta separación no es arbitraria; garantiza que cada nodo rojo tenga únicamente vecinos negros, y que cada nodo negro tenga solo vecinos rojos. Esto permite eliminar dependencias y volver paralelizable un método que originalmente era secuencial.

De esta forma, primero se pueden actualizar todos los nodos rojos en paralelo, y luego todos los nodos negros (o al revés). Gracias a esto, cada conjunto de nodos puede procesarse sin esperar al otro, permitiendo paralelizar Gauss–Seidel de manera efectiva.
```cpp
for (int fase = 0; fase < 2; fase++) {   // 0 = rojos, 1 = negros
    #pragma omp parallel for collapse(2) reduction(max:maxDiff) schedule(static)
    for (int bijy = 0; bijy < nBlocksY; bijy++) {
        for (int bxi = 0; bxi < nBlocksX; bxi++) {
            ...
            for (int j = jStart; j < jEnd; j++) {
                for (int i = iStart; i < iEnd; i++) {
                    if ((i + j) % 2 != fase) continue; // Red-Black

                    double newValue = ...
                    maxDiff = max(maxDiff, fabs(newValue - current[idx]));
                    current[idx] = newValue;
                }
            }
        }
    }
}
```

## Directivas OpenMP

### `#pragma omp parallel for collapse(2)`

Este pragma indica que los bucles externos (los que recorren los bloques en las direcciones X e Y) se deben paralelizar. El uso de `collapse(2)` le dice a OpenMP que combine ambos bucles en uno solo, creando un conjunto de iteraciones más grande y uniforme. Esto ayuda a repartir mejor los bloques entre los hilos y evita que algunos queden sin trabajo cuando una dimensión tiene pocos bloques.

- **Sin collapse**: solo se paraleliza por filas de bloques.
- **Con collapse**: todos los bloques se reparten equitativamente.

Sin collapse, es común que algunos hilos se queden esperando mientras otros aún siguen trabajando.

### `schedule(static)`

El modo `static` se usa porque:

- Todos los bloques requieren prácticamente el mismo esfuerzo computacional.
- El tamaño del trabajo por bloque es muy similar.
- Se evita el costo adicional del planificador `dynamic`, que introduce más sincronizaciones.

Con esto se logra una distribución rápida, simple y eficiente del cálculo entre hilos.

### `reduction(max:maxDiff)`

En cada iteración se debe calcular la norma infinita del error:

$$
\max|u^{(k+1)} - u^{(k)}|
$$

Como varios hilos contribuyen a este cálculo al mismo tiempo, se utiliza una reducción de tipo `max`. Esto permite que cada hilo compute su propio valor local sin interferir con los demás y que al final todos se combinen de forma segura para obtener el máximo global, evitando condiciones de carrera.

## División en Bloques (Tiling)

Además del coloreo Red-Black, la malla también se divide en bloques rectangulares (por ejemplo, de 32 x 32 celdas). Cada bloque se asigna a un hilo mediante OpenMP. Esto aporta varios beneficios importantes:

- Se reducen los accesos a memoria RAM
- Los datos cercanos se reutilizan dentro del bloque
- Los bloques se distribuyen de forma balanceada entre los hilos
- Se evita que varios hilos trabajen sobre datos demasiado próximos, reduciendo el riesgo de *false sharing*

La combinación del tiling con el esquema Red–Black genera un paralelismo muy eficiente que escala bien con el número de núcleos del procesador.
```cpp
const int BX = 32, BY = 32;
const int nBlocksX = (N + BX - 1) / BX;
const int nBlocksY = (N + BY - 1) / BY;

int jStart = bjy * BY, jEnd = min(N, jStart + BY);
int iStart = bxi * BX, iEnd = min(N, iStart + BX);
```

Así, cada hilo trabaja sobre un bloque bien delimitado, cuyos datos permanecen dentro de las caches del procesador durante casi toda la operación, lo que incrementa significativamente el rendimiento general.

