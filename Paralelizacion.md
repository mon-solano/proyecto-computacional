# Paralelización Red–Black con Bloques (Tiling)

En la resolución numérica de la ecuación de calor en dos dimensiones mediante el método de Euler implícito, cada paso de tiempo requiere resolver un sistema lineal grande. Este sistema se resuelve con Gauss–Seidel, pero su versión clásica es totalmente secuencial, ya que cada punto depende de las actualizaciones recientes de sus vecinos. Esto hace que el método sea lento cuando la malla es grande. Para aprovechar los procesadores multinúcleo y reducir los tiempos de cómputo, se combinan dos técnicas muy utilizadas en paralelización: el esquema Red–Black Gauss–Seidel y la división de la malla en bloques o *tiles*.

La idea del método Red–Black es dividir la malla como si fuera un tablero de ajedrez: los nodos donde \(i + j\) es par se marcan como rojos, y los nodos donde \(i + j\) es impar como negros. Esta separación no es arbitraria; garantiza que cada nodo rojo tenga únicamente vecinos negros, y que cada nodo negro tenga solo vecinos rojos. Esto permite eliminar dependencias y volver paralelizable un método que originalmente era secuencial.

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


#pragma omp parallel for collapse(2)

Este pragma indica que los bucles externos (los que recorren los bloques en las direcciones X y Y) se deben paralelizar. El uso de collapse(2) le dice a OpenMP que combine ambos bucles en uno solo, creando un conjunto de iteraciones más grande y uniforme. Esto ayuda a repartir mejor los bloques entre los hilos y evita que algunos queden sin trabajo cuando una dimensión tiene pocos bloques.

    Sin collapse: solo se paraleliza por filas de bloques.

    Con collapse: todos los bloques se reparten equitativamente.

Sin collapse, es común que algunos hilos se queden esperando mientras otros aún siguen trabajando.

schedule(static)

El modo static se usa porque:

    Todos los bloques requieren prácticamente el mismo esfuerzo computacional.

    El tamaño del trabajo por bloque es muy similar.

    Se evita el costo adicional del planificador dynamic, que introduce más sincronizaciones.

Con esto se logra una distribución rápida, simple y eficiente del cálculo entre hilos.

En cada iteración se debe calcular la norma infinita del error:

\[
\max|u^{(k+1)} - u^{(k)}|
\]

Como varios hilos contribuyen a este cálculo al mismo tiempo, se utiliza una reducción de tipo `max`. Esto permite que cada hilo compute su propio valor local sin interferir con los demás y que al final todos se combinen de forma segura para obtener el máximo global, evitando condiciones de carrera.


Además del coloreo Red-Black, la malla también se divide en bloques rectangulares (por ejemplo, de 32 x 32 celdas). Cada bloque se asigna a un hilo mediante OpenMP. Esto aporta varios beneficios importantes: se reducen los accesos a memoria RAM, los datos cercanos se reutilizan dentro del bloque, los bloques se distribuyen de forma balanceada entre los hilos y se evita que varios hilos trabajen sobre datos demasiado próximos, reduciendo el riesgo de false sharing.

La combinación del tiling con el esquema Red–Black genera un paralelismo muy eficiente que escala bien con el número de núcleos del procesador.


```cpp

const int BX = 32, BY = 32;
const int nBlocksX = (N + BX - 1) / BX;
const int nBlocksY = (N + BY - 1) / BY;

int jStart = bjy * BY, jEnd = min(N, jStart + BY);
int iStart = bxi * BX, iEnd = min(N, iStart + BX);

Así, cada hilo trabaja sobre un bloque bien delimitado, cuyos datos permanecen dentro de las caches del procesador durante casi toda la operación, lo que incrementa significativamente el rendimiento general.
