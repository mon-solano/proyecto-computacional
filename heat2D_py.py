import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve
import time
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

class HeatEquation2D:

    def __init__(self, N, M, T, a, b, c, condiciones_frontera):
        self.N = N # numero de puntos interiores en x y y
        self.M = M # numero de pasos temporales
        self.T = T # tiempo final de simulacion
        self.a = a # longitud del dominio en x
        self.b = b # longitud del dominio en y
        self.c = c # coeficiente de difusion termica

        # Paso espacial y temporal, todo esto esta definido en la teoria
        self.hx = a / (N + 1) # N + 1 pq si se tienen N puntos interiores, hay que dividir el dominio en N+1 segmentos
        self.hy = b / (N + 1)
        self.k = T / M

        # Parámetros de discretización, esto se deriva en la teoria, por la estabilidad del metodo, no tiene que ser algun valor especifico
        self.mux = (c**2 * self.k) / self.hx**2
        self.muy = (c**2 * self.k) / self.hy**2

        # Malla espacial sin incluir las fronteras
        self.x = np.linspace(self.hx, a - self.hx, N)
        self.y = np.linspace(self.hy, b - self.hy, N)
        self.X, self.Y = np.meshgrid(self.x, self.y)
        self.condiciones_frontera = condiciones_frontera

        # Construir matriz del sistema
        self.A = self.MatrizLaplaciano()

        print(f"Configuración del solver:")
        print(f"  N = {N} (puntos interiores por dimensión)")
        print(f"  M = {M} (pasos de tiempo)")
        print(f"  hx = {self.hx:.6f} (paso espacial en x)")
        print(f"  hy = {self.hy:.6f} (paso espacial en y)")
        print(f"  k = {self.k:.6f} (paso temporal)")
        print(f"  μx = {self.mux:.6f}")
        print(f"  μy = {self.muy:.6f}")
        print(f"  Tamaño del sistema: {N*N} × {N*N}")
        print(f"\nCondiciones de frontera:")
        for lado, cond in self.condiciones_frontera.items():
            print(f"  {lado}: {cond['type']} = {cond['valor']}")

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
                g_top = self.condiciones_frontera['top']['valor']
                b[idx_top] += self.muy * g_top

            elif self.condiciones_frontera['top']['type'] == 'neumann':
                g_top = self.condiciones_frontera['top']['valor']
                b[idx_top] += 2.0 * self.muy * self.hy * g_top

        return b


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

# condiciones iniciales con temp. inicial de 100 grados celsius

def DistribucionGaussiana(x, y):
    """Gaussiana centrada en (0.5, 0.5)"""
    return 100 * np.exp(-50 * ((x - 0.5)**2 + (y - 0.5)**2))

def Anillo(x, y):
    """Forma de anillo"""
    r = np.sqrt((x - 0.5)**2 + (y - 0.5)**2)
    return 100 * np.exp(-200 * (r - 0.3)**2)

def Diagonal(x, y):
    """Diagonal"""
    return 100 * np.exp(-100 * (x - y)**2)

def FormaU(x, y):
    """Forma de letra U"""
    left = np.exp(-200 * ((x - 0.25)**2 + (y - 0.5)**2))
    right = np.exp(-200 * ((x - 0.75)**2 + (y - 0.5)**2))
    bottom = np.exp(-200 * ((x - 0.5)**2 + (y - 0.25)**2))
    return 100 * (left + right + bottom)

def CuatroEsquinas(x, y):
  # cuatro esquinas calientes
    sigma = 0.08 # anchura de los puntos calientes
    return 100 * (
        np.exp(-((x - 0.1)**2 + (y - 0.1)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.9)**2 + (y - 0.1)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.1)**2 + (y - 0.9)**2) / (2 * sigma**2)) +
        np.exp(-((x - 0.9)**2 + (y - 0.9)**2) / (2 * sigma**2)))

def analyze_min_max(soluciones, tiempos):
    print("\nAnálisis de temperatura:")
    for i, (sol, t) in enumerate(zip(soluciones, tiempos)):
        u_max = sol.max()
        u_min = sol.min()
        print(f"  t={t:.4f} | Temp. máx: {u_max:.4f} | Temp. mín: {u_min:.4f}")

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


def crear_animacion(soluciones, tiempos, titulo, condiciones_frontera, intervalo=500):
    """
    Crea una animación de la evolución de la ecuación de calor.

    Parámetros:
    -----------
    soluciones : list of ndarray
        Lista con las soluciones en cada tiempo guardado
    tiempos : list of float
        Lista con los tiempos correspondientes
    titulo : str
        Título de la animación
    boundary_conditions : dict
        Condiciones de frontera usadas
    intervalo : int
        Milisegundos entre frames (más alto = más lento)
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    vmin = min(sol.min() for sol in soluciones)
    vmax = max(sol.max() for sol in soluciones)

    im = ax.imshow(soluciones[0], cmap='hot', interpolation='bilinear',
                   origin='lower', vmin=vmin, vmax=vmax, extent=[0, 1, 0, 1])

    cbar = plt.colorbar(im, ax=ax, label='Temperatura')

    # Texto con información dinámica
    texto_info = ax.text(0.02, 0.98, '', transform=ax.transAxes,
                        verticalalignment='top', bbox=dict(boxstyle='round',
                        facecolor='white', alpha=0.8), fontsize=10)

    # Información de condiciones de frontera (estática)
    bc_text = "Condiciones de frontera:\n"
    for lado, cond in condiciones_frontera.items():
        bc_text += f"{lado}: {cond['type']} = {cond['valor']}\n"

    texto_bc = ax.text(0.98, 0.98, bc_text, transform=ax.transAxes,
                      verticalalignment='top', horizontalalignment='right',
                      bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7),
                      fontsize=8)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title(titulo, fontsize=14, fontweight='bold', pad=20)

    def actualizar(frame):
        im.set_array(soluciones[frame])

        t = tiempos[frame]
        temp_max = soluciones[frame].max()
        temp_min = soluciones[frame].min()

        info = f't = {t:.4f}\nMax = {temp_max:.2f}\nMin = {temp_min:.2f}'
        texto_info.set_text(info)

        return [im, texto_info]

    anim = FuncAnimation(fig, actualizar, frames=len(soluciones),
                        interval=intervalo, blit=True, repeat=True)

    plt.tight_layout()
    plt.close()  # Cerrar para evitar mostrar frame estático
    return anim

if __name__ == "__main__":
    print("=" * 70)
    print("ECUACIÓN DE CALOR EN 2D")
    print("=" * 70)

    # Seleccionar condición inicial
    condicion_inicial = seleccionar_condicion_inicial()

    # Solicitar condiciones de frontera al usuario
    condicion_frontera = solicitar_condiciones_frontera()

    # Resolver con las condiciones especificadas
    print("Iniciando simulación...")

    solver = HeatEquation2D(N=250, M=500, T=0.65, a=1, b=1, c=1,
                            condiciones_frontera=condicion_frontera)

    soluciones, tiempos = solver.solve(condicion_inicial, save_every=5)

    analyze_min_max(soluciones, tiempos)

    print("Generando animación...")

    # crear animación
    nombre_condicion = {
        DistribucionGaussiana: "Distribución Gaussiana",
        Anillo: "Anillo",
        Diagonal: "Diagonal",
        FormaU: "Forma de U",
        CuatroEsquinas: "Cuatro esquinas",
    }.get(condicion_inicial, "Condición inicial personalizada")

    anim = crear_animacion(soluciones, tiempos,
                          f"Ecuación de Calor 2D - {nombre_condicion}",
                          condicion_frontera, intervalo=500)

    print("\nAnimación generada. Mostrando...")

    # esto es para mostrar como animacion dentro de colab
    display(HTML(anim.to_jshtml()))

    print("Simulación completada.")
