Este método utiliza condiciones de frontera de Dirichlet no homogéneas para poder resolver la ecuación de calor en dos dimensiones. Aquí se explica el procedimiento con sus diferentes pasos a modo resumen, cómo validarlo usando MATLAB y una idea de cómo aplicarlo a nuestro proyecto.

**1 - ESTRUCTURA DE LA ECUACIÓN DE CALOR**

Lo que queremos es encontrar la temperatura $u(x,y,t)$ en una región cuadrada del plano conforme va pasando el tiempo. Esta temperatura va a seguir la siguiente ecuación:

$$
\frac{\partial u}{\partial t} = \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} + f(t,x,y)
$$

donde:
- $\frac{\partial u}{\partial t}$ es el cambio de temperatura con respecto al tiempo.
- $\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}$ representa cómo se esparce la temperatura en el espacio.
- $f(t,x,y)$ es una fuente de calor **(en nuestro caso es cero porque no hay fuente externa)**.

Además de esto hay condiciones **iniciales** (cómo empieza la temperatura en todo el dominio) y de **frontera** (qué temperatura hay en los bordes del dominio) que siguen la sisguiente estructura:

$$
\left\{
\begin{array}{ll}
u(0,x,y) = u_0(x,y), & \text{(condición inicial)} \\[6pt]
u(t,x,y) = g(t,x,y), & \text{(condición de frontera)}
\end{array}
\right.
$$

Ahora bien, la implementación del MDF consiste en reemplazar las derivadas de la ecuación por aproximaciones usando puntos en una malla. Para esto se realizan una serie de pasos explicados en las siguientes secciones.

**2 - DISCRETIZACIÓN DEL DOMINIO**

Se considera un dominio cuadrado unitario $Q = [0,1] \times [0,1]$ y se divide una malla en puntos discretizando el tiempo y el espacio:
- **Tiempo:** $M + 1$ instantes.
- **Espacio:** $N + 2$ puntos en cada dirección (incluyendo bordes).

Esto generará una malla cartesiana de $(N+2)^2$ puntos.

Ahora se definen los pasos de las variables:
- $h = \frac{1}{N + 1}$: paso en el espacio.
- $k = \frac{T}{M}$: paso en el tiempo.

**3 - APROXIMACIÓN DE LAS DERIVADAS**

Primero, se usa un método llamado ´´esquema de 5 puntos´´ para poder aproximar las **derivadas espaciales**, o bien, el Laplaciano:

$$
\Delta u(t_n, x_i, y_j) \approx \frac{1}{h^2}
\left( u_{i+1,j} + u_{i-1,j} + u_{i,j+1} + u_{i,j-1} - 4u_{i,j} \right)
$$

Esto utiliza los puntos en las 4 direcciones cardinales del punto $(i,j)$.

Para la derivada temporal se usa el método de Euler implícito (backward Euler):

$$
\frac{u_{i,j}^n - u_{i,j}^{n-1}}{k} \approx derivada\space temporal
$$

Este método es incondicionalmente estable, lo que significa que no importa el tamaño de h y k, el método no explota.

**4 - MONTAJE DEL SISTEMA LINEAL**

Al combinar las aproximaciones, obtenemos una ecuación para cada punto interior de la malla:

$$
- u_{i-1,j}^n - u_{i,j-1}^n + (4 + \mu) u_{i,j}^n - u_{i+1,j}^n - u_{i,j+1}^n = h^2 f_{i,j}^n + \mu u_{i,j}^{n-1}
$$

donde $\mu = \dfrac{h^2}{k}$

Para cada paso de tiempo que se realice, se resuelve un sistema lineal de la siguiente forma:

$$
Au^n = b^n
$$

$A$ es una matriz constante de tamaño $N^2 \times N^2$ independiente del tiempo y construida con productos de Kronecker:

$$
A = I \otimes T + S \otimes I
$$

donde:
- $T$ es tridiagonal con $(4+\mu)$ en la diagonal $(y-1)$ en las sub/superdiagonales.
- $S$ es tridiagonal con $-1$ en las sub/superdiagonales.
- $I$ es la identidad de orden $N$.

En MATLAB se puede utilizar el siguiente código para poder contruir la matriz $A$:

```matlab
function [A] = Laplacian2DMatrix(N, mu)
  I = speye(N);
  e = -1 * ones(N, 1);
  T = spdiags([e, -(4 + mu) * e, e], [-1, 0, 1], N, N);
  S = spdiags([e, e], [-1, 1], N, N);
  A = kron(I, T) + kron(S, I);
end
```

El vector $u^n$ es un vector de incógnitas que representa la temperatura en todos los puntos interiores en el tiempo $t_n$.

Por último, el vector $b^n$ en el lado derecho de la ecuación se construye para cada tiempo $t_n$ usando:
- El término fuente $f(t_n, x_i, y_j)$
- Las condiciones de frontera $g(t_n, x, y)$
- La solución del paso anterior $u^{n-1}$

Este se puede hacer con una rutina que recorra todos los puntos y sume las contribuciones, como se aplica en el siguiente código:

```matlab
function [b] = Lapacian2DRhs(N, tn, x, y, h, f, g, mu, u)
  b = zeros(N^2, 1);
  k = N^2 - N;
  v = y(1); w = y(N);
  for i = 1:N
    z = x(i);
    b(i) = h^2 * f(tn, z, v) + g(tn, z, 0);
    b(k + i) = h^2 * f(tn, z, w) + g(tn, z, 1);
  end
  b(1) = b(1) + g(tn, 0, v);
  b(N) = b(N) + g(tn, 1, v);
  b(k+1) = b(k+1) + g(tn, 0, w);
  b(k+N) = b(k+N) + g(tn, 1, w);
  k = N;
  for j = 2:N-1
    v = y(j);
    for i = 1:N
      b(k + i) = h^2 * f(tn, x(i), v);
    end
    b(k+1) = b(k+1) + g(tn, 0, v);
    k = k + N;
    b(k) = b(k) + g(tn, 1, v);
  end
  b = b + mu * u;
end
```

**5 - EVOLUCIÓN TEMPORAL DE LA SOLUCIÓN**

Se inicializa la solución con la condición inicial $u_0(x, y)$ y se avanza en el tiempo con un bucle desde $n=1$ hasta $M$ siguiendo el siguiente orden:
- Se calcula $b^n$
- Se resuelve $u^n = A^{-1} b^n$
- Se guarda $u^n$ para el siguiente paso

Al final lo que se obtiene es la solución aproximada en todos los tiempos.

Esto se puede lograr con el siguiente código:

```matlab
function [u] = HeatEquation2DScheme(T, N, M, f, g, u0)
  h = 1 / (N + 1);
  x = h : h : (1 - h);
  y = x;
  k = T / M;
  t = k : k : T;
  u = zeros(N^2, 1);
  for j = 1:N
    for i = 1:N
      u((j-1)*N + i) = u0(x(i), y(j));
    end
  end
  mu = h^2 / k;
  A = Laplacian2DMatrix(N, mu);
  for n = 1:M
    b = Lapacian2DRhs(N, t(n), x, y, h, f, g, mu, u);
    u = A \ b;
  end
end
```

**6 - VALIDACIÓN DEL MÉTODO**

Para validar el código anterior se propone una solución exacta conocida:

$$
u(t,x,y)=sin(t+πx)⋅sin(t+πy)
$$

Y se calcula la fuente $f(t, x, y)$ que hace que esta función cumpla la ecuación de calor.
Luego, se compara la solución numérica con la exacta para verificar la precisión del método.

**En nuestro caso no existe ninguna fuente externa de energía, por lo que si se llegara a realizar la comparación, la única forma de que ambas soluciones estén correctas es que ambas den cero.**

##¿CÓMO APLICAR ESTE MÉTODO A NUESTRO PROYECTO?

Nuestra ecuación se define como:

$$
\frac{\partial u}{\partial t} = c^2 (\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2})
$$

Es básicamente lo mismo que se explica anteriormente pero con $f=0$ y un coeficiente de difusión térmica $c^2$.

Para poder aplicar el método ocupamos hacer lo siguiente:
- Redefinir el dominio como $x\in [0,a]$ y $y\in [0,b]$
- Definir los pasos en el espacio como:
$$
h_x = \frac{a}{N + 1}\space h_y = \frac{b}{N + 1}
$$
Si tomamos un $a=b$ entonces se puede usar un mismo h.
- Si $h_x\neq h_y$ entonces se tiene que modificar el Laplaciano para que quede de la siguiente forma:
$$
\Delta u\approx \frac{u_{i+1,j}-2u_{i,j}+u_{i-1,j}}{h_x^2}+\frac{u_{i,j+1}-2u_{i,j}+u_{i,j-1}}{h_y^2}
$$
- Se debe incluir el coeficiente $c^2$ multiplicándolo por todo el Laplaciano. Hay que tomar en cuenta que esto afecta a la matriz $A$ ya que cambia los valores de $\mu$ y los coeficientes.
- Hay que definir una condición de frontera $g(t,x,y)$ según el caso (puede ser constante, cero o una función).
- Se define una condición inicial tomando a $u_0(x,y)$ como la temperatura inicial.
- Se toma que $\mu = \frac{h^2}{kc^2}$
