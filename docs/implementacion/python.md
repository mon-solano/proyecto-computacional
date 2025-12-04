# Implementación en Python

La implementación en Python utiliza NumPy para operaciones matriciales y SciPy para resolver el sistema lineal disperso de manera eficiente.

## Arquitectura del Código

### Clase Principal: `HeatEquation2D`

```python
class HeatEquation2D:
    def __init__(self, N=150, M=200, T=0.5, a=1.0, b=1.0, c=1.0,
                 condiciones_frontera=None)
    def MatrizLaplaciano(self)
    def VectorRHS(self, tn, u_prev)
    def obtener_valor_frontera(self, lado, t, pos)
    def solve(self, u0_func, save_every=1)
```

## Componentes Principales

### 1. Inicialización
Aqui simplemente se definen los parametros con los que se van a trabajar y se definen como atributos de la clase:

```python
def __init__(self, N=150, M=200, T=0.5, a=1.0, b=1.0, c=1.0,
             condiciones_frontera=None):
    self.N = N
    self.M = M
    self.T = T
    self.a = a
    self.b = b
    self.c = c
    
    # Pasos espaciales y temporales
    self.hx = a / (N + 1)
    self.hy = b / (N + 1)
    self.k = T / M
    
    # Parámetros de discretización
    self.mux = (c**2 * self.k) / self.hx**2
    self.muy = (c**2 * self.k) / self.hy**2
    
    # Malla espacial (solo puntos interiores)
    self.x = np.linspace(self.hx, a - self.hx, N)
    self.y = np.linspace(self.hy, b - self.hy, N)
    self.X, self.Y = np.meshgrid(self.x, self.y)
    
    # Construir matriz del sistema
    self.A = self.MatrizLaplaciano()
```


### 2. Construcción de la Matriz del Sistema
Se elige construir la matriz A en formato LIL ya que es un formato eficiente para construir matrices dispersas elemento por elemento. El elemento idx proyecta los elementos \((i,j)\) en una dimension ordenandolos por fila. Ademas, se guarda como csr para optimizar las operaciones algebraicas.

```python
def MatrizLaplaciano(self):
        N = self.N
        A = lil_matrix((N * N, N * N))

        # Pre-calcular tipos de condiciones para cada borde
        left_type = self.condiciones_frontera['left']['type']
        right_type = self.condiciones_frontera['right']['type']
        bottom_type = self.condiciones_frontera['bottom']['type']
        top_type = self.condiciones_frontera['top']['type']

        for j in range(N):
            for i in range(N):
                idx = j * N + i

                # Diagonal: siempre 1 + 2μx + 2μy
                A[idx, idx] = 1.0 + 2.0 * self.mux + 2.0 * self.muy

                # Vecinos en x
                # Está en borde izquierdo
                if i == 0:
                    if left_type == 'neumann':
                        # En borde izquierdo con Neumann
                        # No hay coeficiente para vecino izquierdo (fantasma)
                        # Pero el vecino derecho (si existe) tiene -2μx
                        if i < N - 1:  # Existe vecino derecho
                            A[idx, idx + 1] = -2.0 * self.mux
                    else:
                        # Dirichlet en borde izquierdo
                        # Coeficiente normal para vecino derecho
                        if i < N - 1:
                            A[idx, idx + 1] = -self.mux
                # Está en borde derecho
                elif i == N - 1:
                    if right_type == 'neumann':
                        # En borde derecho con Neumann
                        # Vecino izquierdo tiene -2μx
                        A[idx, idx - 1] = -2.0 * self.mux
                    else:
                        # Dirichlet en borde derecho
                        # Coeficiente normal para vecino izquierdo
                        A[idx, idx - 1] = -self.mux
                else:
                    # Punto interior en x
                    A[idx, idx - 1] = -self.mux  # Vecino izquierdo
                    A[idx, idx + 1] = -self.mux  # Vecino derecho

                # Vecinos en Y

                # Está en borde inferior
                if j == 0:
                    if bottom_type == 'neumann':
                        # En borde inferior con Neumann
                        # Vecino superior (si existe) tiene -2μy
                        if j < N - 1:
                            A[idx, idx + N] = -2.0 * self.muy
                    else:
                        # Dirichlet en borde inferior
                        if j < N - 1:
                            A[idx, idx + N] = -self.muy
                # Está en borde superior
                elif j == N - 1:
                    if top_type == 'neumann':
                        # En borde superior con Neumann
                        # Vecino inferior tiene -2μy
                        A[idx, idx - N] = -2.0 * self.muy
                    else:
                        # Dirichlet en borde superior
                        A[idx, idx - N] = -self.muy
                else:
                    # Punto interior en y
                    A[idx, idx - N] = -self.muy  # Vecino inferior
                    A[idx, idx + N] = -self.muy  # Vecino superior

        return A.tocsr()
```

### 3. Construcción del Vector RHS
Aqui simplemente se construye el vector segun el tipo de condicion y el borde que se tenga.

```python
def VectorRHS(self, tn, u_prev):
        """
        aqui se construye el vector del lado derecho

        se dedujo que:
        - borde izquierdo Neumann:  b = uⁿ - 2μx·hx·g
        - borde derecho Neumann:    b = uⁿ + 2μx·hx·g
        - borde inferior Neumann:   b = uⁿ - 2μy·hy·g
        - borde superior Neumann:   b = uⁿ + 2μy·hy·g
        - Dirichlet:                b = uⁿ + μ·g

        esto es lo que se desarrolla aca
        """
        N = self.N
        b = u_prev.copy()  # empezamos con las temperaturas anteriores (o sea, uⁿ)

        # bordes izquierdo y derecho (i fijo)
        for j in range(N):
            y_pos = self.y[j]  # posicion en y

            # borde izquierdo (i=0)
            idx_left = j * N  # j*N + 0
            if self.condiciones_frontera['left']['type'] == 'dirichlet':
                # Dirichlet: u = g en el borde
                g_left = self.condiciones_frontera['left']['valor']
                # se añade el término necesario
                b[idx_left] += self.mux * g_left

            elif self.condiciones_frontera['left']['type'] == 'neumann':
                # Neumann: ∂u/∂x = g en el borde
                g_left = self.condiciones_frontera['left']['valor']
                # se añade el término necesario
                b[idx_left] -= 2.0 * self.mux * self.hx * g_left

            # borde derecho (i=N-1), todo es la misma lógica que para el borde
            # izq
            idx_right = j * N + (N - 1)
            if self.condiciones_frontera['right']['type'] == 'dirichlet':
                g_right = self.condiciones_frontera['right']['valor']
                b[idx_right] += self.mux * g_right

            elif self.condiciones_frontera['right']['type'] == 'neumann':
                g_right = self.condiciones_frontera['right']['valor']
                b[idx_right] += 2.0 * self.mux * self.hx * g_right

        # bordes inferior y superior (j fijo), igual todo es con la misma
        # lógica y se usan las mismas fórmulas
        for i in range(N):
            x_pos = self.x[i]  # posición en x

            # borde inferior (j=0)
            idx_bottom = i  # 0*N + i = i
            if self.condiciones_frontera['bottom']['type'] == 'dirichlet':
                g_bottom = self.condiciones_frontera['bottom']['valor']
                b[idx_bottom] += self.muy * g_bottom

            elif self.condiciones_frontera['bottom']['type'] == 'neumann':
                g_bottom = self.condiciones_frontera['bottom']['valor']
                b[idx_bottom] -= 2.0 * self.muy * self.hy * g_bottom

            # borde superior (j=N-1)
            idx_top = (N - 1) * N + i
            if self.condiciones_frontera['top']['type'] == 'dirichlet':
                g_top = self.obtener_valor_frontera['top']['valor']
                b[idx_top] += self.muy * g_top

            elif self.condiciones_frontera['top']['type'] == 'neumann':
                g_top = self.condiciones_frontera['top']['valor']
                b[idx_top] += 2.0 * self.muy * self.hy * g_top

        return b
```


### 4. Solver Principal

Se inicializa el vector \(u0\) con la condicion inicial evaluada en cada puntos y se almacena en la lista de soluciones.

Despues hay un bucle temporal en el que para cado paso n se calcula el tiempo actual, ademas construye el vector b basandose en \(u_prev\) y en las condiciones de frontera en el tiempo actual.

Con esto ya armado, se resuelve para u_actual con la funcion spsolve, que es una funcion de SciPy para resolver sistemas lineales dispersos.

```python
def solve(self, u0_func, save_every=1):
        # resuelve la ecuación de calor
        N = self.N

       # condicion inicial
        u0 = np.zeros(N * N)
        for j in range(N):
            for i in range(N):
                idx = j * N + i
                u0[idx] = u0_func(self.x[i], self.y[j])

        soluciones = [u0.reshape(N, N).copy()]
        tiempos = [0.0]
        u_prev = u0

        print(f"Resolviendo...")

        # evolución temporal
        for n in range(1, self.M + 1):
            tn = n * self.k  # tiempo actual

            # se construye el vector del lado derecho
            b = self.VectorRHS(tn, u_prev)

            # se resuelve sistema lineal A·u^{n+1} = b
            u_actual = spsolve(self.A, b)

            if n % save_every == 0 or n == self.M:
                soluciones.append(u_actual.reshape(N, N).copy())
                tiempos.append(tn)

            u_prev = u_actual

        return soluciones, tiempos
```

## Condiciones 

Se ponen cinco condiciones iniciales con una temperatura inicial de 100 grados celsius:

### 1. Distribución Gaussiana

```python
def DistribucionGaussiana(x, y):
    """Pico gaussiano centrado en (0.5, 0.5)"""
    return 100 * np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))
```

### 2. Anillo

```python
def Anillo(x, y):
    """Distribución en forma de anillo"""
    r = np.sqrt((x - 0.5)**2 + (y - 0.5)**2)
    return 100 * np.exp(-200 * (r - 0.3)**2)
```

### 3. Diagonal

```python
def Diagonal(x, y):
    """Franja diagonal caliente"""
    return 100 * np.exp(-100 * (x - y)**2)
```

### 4. Forma de U

```python
def FormaU(x, y):
    """Tres puntos calientes formando U"""
    left = np.exp(-200 * ((x - 0.25)**2 + (y - 0.5)**2))
    right = np.exp(-200 * ((x - 0.75)**2 + (y - 0.5)**2))
    bottom = np.exp(-200 * ((x - 0.5)**2 + (y - 0.25)**2))
    return 100 * (left + right + bottom)
```

### 5. Cuatro Esquinas

```python
def CuatroEsquinas(x, y):
    """Cuatro puntos calientes en las esquinas"""
    sigma = 0.08
    return 100 * (
        np.exp(-((x - 0.1)**2 + (y - 0.1)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.9)**2 + (y - 0.1)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.1)**2 + (y - 0.9)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.9)**2 + (y - 0.9)**2) / (2 * sigma**2))
    )
```

Ademas de esto, se implementa una funcion que imprime las temperaturas maximas y minimas en el sistema en un determinado tiempo.

```python
def analyze_min_max(soluciones, tiempos):
    print("\nAnálisis de temperatura:")
    for i, (sol, t) in enumerate(zip(soluciones, tiempos)):
        u_max = sol.max()
        u_min = sol.min()
        print(f"  t={t:.4f} | Temp. máx: {u_max:.4f} | Temp. mín: {u_min:.4f}")
```

## Menu interactivo 

Se implementaron funciones que permiten al usuario elegir las condiciones iniciales y las condiciones de frontera para cada uno de los bordes de la malla.

```python
def seleccionar_condicion_inicial():
    """
    Permite al usuario seleccionar una condición inicial.
    """
    print("Selección de condición inicial")
    print("  1. Distribución Gaussiana")
    print("  2. Anillo")
    print("  3. Diagonal")
    print("  4. Forma de U")
    print("  5. Cuatro esquinas")

    condiciones = {
        1: ("Distribución Gaussiana", DistribucionGaussiana),
        2: ("Anillo", Anillo),
        3: ("Diagonal", Diagonal),
        4: ("Forma de U", FormaU),
        5: ("Cuatro esquinas", CuatroEsquinas)
    }

    while True:


      while True:
        try:
            opcion = int(input("\nSeleccione una opción (1-5): ").strip())
            if opcion in condiciones:
                nombre, funcion = condiciones[opcion]
                print(f"\nSeleccionado: {nombre}")
                return funcion
            else:
                print("Error: Ingrese un número entre 1 y 4")
        except ValueError:
            print("Error: Ingrese un número válido")


def solicitar_condiciones_frontera():
    """
    Solicita al usuario las condiciones de frontera para cada lado del dominio.
    """
    print("Configuración de condiciones de frontera")
    print("\nTipos de condición:")
    print("  Dirichlet: temperatura fija en la frontera")
    print("  Neumann: flujo de calor en la frontera (0 = aislado)")

    condiciones_frontera = {}
    lados = ['left', 'right', 'bottom', 'top']
    lado_nombres = {
        'left': 'Izquierdo',
        'right': 'Derecho',
        'bottom': 'Inferior',
        'top': 'Superior'
    }

    for lado in lados:
        print(f"\n--- Borde {lado_nombres[lado]} ---")

        while True:
            tipo = input(f"Tipo de condición (dirichlet/neumann): ").strip().lower()
            if tipo in ['dirichlet', 'neumann', 'd', 'n']:
                if tipo == 'd':
                    tipo = 'dirichlet'
                elif tipo == 'n':
                    tipo = 'neumann'
                break
            print("Error: Ingrese 'dirichlet' (o 'd') o 'neumann' (o 'n')")

        while True:
            try:
                valor_str = input(f"Valor: ").strip()
                valor = float(valor_str)
                break
            except ValueError:
                print("Error: Ingrese un número válido")

        condiciones_frontera[lado] = {'type': tipo, 'valor': valor}

    return condiciones_frontera
```

## Visualización y Animación

### Función de Animación

```python
def crear_animacion(soluciones, tiempos, titulo, condiciones_frontera, 
                   intervalo=500):
    """
    Crea animación de la evolución de temperatura.
    
    Parámetros:
    -----------
    soluciones : list
        Lista de arrays 2D con soluciones
    tiempos : list
        Tiempos correspondientes
    titulo : str
        Título de la animación
    condiciones_frontera : dict
        Condiciones de frontera usadas
    intervalo : int
        Milisegundos entre frames
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Rango de colores consistente
    vmin = min(sol.min() for sol in soluciones)
    vmax = max(sol.max() for sol in soluciones)
    
    # Imagen inicial
    im = ax.imshow(soluciones[0], cmap='hot', interpolation='bilinear',
                   origin='lower', vmin=vmin, vmax=vmax, 
                   extent=[0, 1, 0, 1])
    
    # Barra de color
    cbar = plt.colorbar(im, ax=ax, label='Temperatura')
    
    # Función de actualización
    def actualizar(frame):
        im.set_array(soluciones[frame])
        # Actualizar información de tiempo y temperatura
        # ...
        return [im]
    
    anim = FuncAnimation(fig, actualizar, frames=len(soluciones),
                        interval=intervalo, blit=True, repeat=True)
    
    return anim
```

---

**Próximo paso**: Revisa [Ejemplos de Uso](../ejemplos/condiciones-iniciales.md) para casos prácticos completos.