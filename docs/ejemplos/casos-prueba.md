# Casos de Prueba

Esta sección presenta casos de prueba completos que combinan diferentes condiciones iniciales y de frontera para ilustrar diversos comportamientos físicos.

## Caso 1: Enfriamiento Total

### Descripción

Una región inicialmente caliente que se enfría completamente debido a fronteras frías.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, DistribucionGaussiana

# Todas las fronteras a temperatura 0
condiciones = {
    'left':   {'type': 'dirichlet', 'valor': 0.0},
    'right':  {'type': 'dirichlet', 'valor': 0.0},
    'bottom': {'type': 'dirichlet', 'valor': 0.0},
    'top':    {'type': 'dirichlet', 'valor': 0.0}
}

solver = HeatEquation2D(
    N=150,
    M=300,
    T=1.0,
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(DistribucionGaussiana, save_every=10)
```

### Resultados Esperados

| Tiempo | Temperatura Máxima | Observación |
|--------|-------------------|-------------|
| t=0.0 | 100.0°C | Pico inicial |
| t=0.2 | ~60°C | Disminución rápida |
| t=0.5 | ~20°C | Cerca de las fronteras |
| t=1.0 | ~5°C | Casi completamente frío |

### Análisis

- El calor **fluye hacia las fronteras** y se disipa
- La temperatura **máxima decrece exponencialmente**
- En \(t \to \infty\), toda la región alcanza 0°C

!!! note "Conservación de energía"
    La energía total del sistema **disminuye** porque se transfiere al exterior a través de las fronteras.

---

## Caso 2: Redistribución con Aislamiento

### Descripción

Una región con distribución no uniforme de calor que se redistribuye sin pérdida de energía total.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, CuatroEsquinas

# Todas las fronteras aisladas
condiciones = {
    'left':   {'type': 'neumann', 'valor': 0.0},
    'right':  {'type': 'neumann', 'valor': 0.0},
    'bottom': {'type': 'neumann', 'valor': 0.0},
    'top':    {'type': 'neumann', 'valor': 0.0}
}

solver = HeatEquation2D(
    N=150,
    M=400,
    T=1.5,
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(CuatroEsquinas, save_every=10)
```

### Resultados Esperados

| Tiempo | Temp. Máxima | Temp. Mínima | Temp. Promedio |
|--------|--------------|--------------|----------------|
| t=0.0 | 100°C | 0°C | ~25°C |
| t=0.3 | ~75°C | ~10°C | ~25°C |
| t=0.8 | ~35°C | ~20°C | ~25°C |
| t=1.5 | ~27°C | ~23°C | ~25°C |

### Análisis

- El calor se **redistribuye uniformemente**
- **No hay pérdida de energía** (fronteras aisladas)
- La temperatura promedio **se conserva** exactamente
- En \(t \to \infty\), la temperatura es **uniforme** en todo el dominio

### Verificación de Conservación

```python
# Calcular energía total en cada tiempo
h = solver.hx * solver.hy
energias = [np.sum(sol) * h for sol in soluciones]

print("Energía total en cada tiempo:")
for t, E in zip(tiempos, energias):
    print(f"  t={t:.3f}: E={E:.6f}")

print(f"\nVariación máxima: {np.ptp(energias):.8f}")
```

**Resultado esperado**: Variación \(< 10^{-6}\)

---

## Caso 3: Flujo Constante por una Frontera

### Descripción

Una región donde entra calor constantemente por una frontera mientras otras permanecen frías.

### Configuración

```python
from heat_equation_2d import HeatEquation2D

def InicialCero(x, y):
    """Región inicialmente fría"""
    return 0.0

# Flujo de calor por la izquierda, otras fronteras frías
condiciones = {
    'left':   {'type': 'neumann', 'valor': 10.0},  # Flujo entrante
    'right':  {'type': 'dirichlet', 'valor': 0.0},
    'bottom': {'type': 'dirichlet', 'valor': 0.0},
    'top':    {'type': 'dirichlet', 'valor': 0.0}
}

solver = HeatEquation2D(
    N=150,
    M=500,
    T=2.0,
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(InicialCero, save_every=10)
```

### Resultados Esperados

- **Inicialmente**: Región completamente fría
- **t ≈ 0.2**: Calentamiento cerca de la frontera izquierda
- **t ≈ 0.8**: Gradiente de temperatura establecido
- **t → ∞**: **Estado estacionario** con gradiente lineal

### Estado Estacionario

En equilibrio (\(\frac{\partial u}{\partial t} = 0\)), la solución satisface:

\[
\nabla^2 u = 0 \quad \text{(Ecuación de Laplace)}
\]

Con las condiciones dadas, se espera un **gradiente aproximadamente lineal** desde la frontera izquierda hacia las otras fronteras.

---

## Caso 4: Condiciones Mixtas Asimétricas

### Descripción

Combinación de diferentes tipos de condiciones en cada frontera para crear patrones de flujo complejos.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, Anillo

# Frontera izquierda caliente, derecha fría, superior/inferior aisladas
condiciones = {
    'left':   {'type': 'dirichlet', 'valor': 80.0},   # Caliente
    'right':  {'type': 'dirichlet', 'valor': 20.0},   # Fría
    'bottom': {'type': 'neumann', 'valor': 0.0},      # Aislada
    'top':    {'type': 'neumann', 'valor': 0.0}       # Aislada
}

solver = HeatEquation2D(
    N=150,
    M=400,
    T=1.0,
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(Anillo, save_every=10)
```

### Resultados Esperados

- **Flujo dominante**: De izquierda (caliente) a derecha (fría)
- **Fronteras superior/inferior**: No afectan el flujo (aisladas)
- **Patrón final**: Gradiente aproximadamente horizontal

### Análisis del Flujo

```python
# Calcular flujo de calor en dirección x
def calcular_flujo_x(sol, hx):
    """Calcula ∂u/∂x usando diferencias centrales"""
    flujo = np.zeros_like(sol)
    flujo[:, 1:-1] = (sol[:, 2:] - sol[:, :-2]) / (2 * hx)
    return flujo

# Para última solución
flujo_x = calcular_flujo_x(soluciones[-1], solver.hx)
print(f"Flujo promedio en x: {np.mean(flujo_x):.4f}")
```

---

## Caso 5: Calentamiento Periódico

### Descripción

Una frontera con temperatura que oscila periódicamente en el tiempo.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, DistribucionGaussiana
import numpy as np

# Frontera izquierda con temperatura oscilante
def temp_oscilante(t, pos):
    """Temperatura que oscila sinusoidalmente"""
    T_media = 50.0
    amplitud = 30.0
    periodo = 0.5
    return T_media + amplitud * np.sin(2 * np.pi * t / periodo)

condiciones = {
    'left':   {'type': 'dirichlet', 'valor': temp_oscilante},
    'right':  {'type': 'dirichlet', 'valor': 0.0},
    'bottom': {'type': 'neumann', 'valor': 0.0},
    'top':    {'type': 'neumann', 'valor': 0.0}
}

solver = HeatEquation2D(
    N=150,
    M=600,
    T=3.0,  # Simular varios periodos
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(DistribucionGaussiana, save_every=10)
```

### Resultados Esperados

- **Onda térmica** propagándose desde la frontera izquierda
- **Atenuación**: La amplitud de la oscilación disminuye con la distancia
- **Desfase**: La oscilación se retrasa a medida que nos alejamos de la frontera

### Análisis de Temperatura en Puntos Específicos

```python
# Extraer temperatura en x = 0.25 (cerca de frontera) y x = 0.75 (lejos)
idx_cerca = int(0.25 * solver.N)
idx_lejos = int(0.75 * solver.N)
j_medio = solver.N // 2

temp_cerca = [sol[j_medio, idx_cerca] for sol in soluciones]
temp_lejos = [sol[j_medio, idx_lejos] for sol in soluciones]

plt.figure(figsize=(10, 6))
plt.plot(tiempos, temp_cerca, label='x=0.25 (cerca)', linewidth=2)
plt.plot(tiempos, temp_lejos, label='x=0.75 (lejos)', linewidth=2)
plt.xlabel('Tiempo')
plt.ylabel('Temperatura')
plt.title('Respuesta a calentamiento periódico')
plt.legend()
plt.grid(True)
plt.show()
```

---

## Caso 6: Test de Convergencia

### Descripción

Verificar que la solución **converge** al refinar la malla espacial y temporal.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, DistribucionGaussiana

condiciones = {
    'left':   {'type': 'dirichlet', 'valor': 0.0},
    'right':  {'type': 'dirichlet', 'valor': 0.0},
    'bottom': {'type': 'dirichlet', 'valor': 0.0},
    'top':    {'type': 'dirichlet', 'valor': 0.0}
}

# Resolver con diferentes resoluciones
resoluciones = [
    (50, 100),    # Grueso
    (100, 200),   # Medio
    (150, 300),   # Fino
    (200, 400)    # Muy fino
]

soluciones_finales = []

for N, M in resoluciones:
    solver = HeatEquation2D(
        N=N, M=M, T=0.5,
        condiciones_frontera=condiciones
    )
    sols, _ = solver.solve(DistribucionGaussiana)
    soluciones_finales.append(sols[-1])
    print(f"N={N}, M={M}: Temp. máx. final = {sols[-1].max():.6f}")
```

### Análisis de Convergencia

```python
# Calcular diferencia entre soluciones sucesivas
for i in range(1, len(soluciones_finales)):
    # Interpolar solución gruesa a malla fina para comparación
    sol_gruesa = soluciones_finales[i-1]
    sol_fina = soluciones_finales[i]
    
    # Diferencia relativa (simplificada)
    diff = abs(sol_fina.max() - sol_gruesa.max()) / sol_fina.max()
    print(f"Diferencia {i-1}→{i}: {diff:.6f} ({diff*100:.3f}%)")
```

**Resultado esperado**: Las diferencias deben **disminuir** al refinar la malla, confirmando convergencia.

---

## Caso 7: Simetría

### Descripción

Verificar que condiciones simétricas producen soluciones simétricas.

### Configuración

```python
from heat_equation_2d import HeatEquation2D, DistribucionGaussiana

# Condiciones completamente simétricas
condiciones = {
    'left':   {'type': 'neumann', 'valor': 0.0},
    'right':  {'type': 'neumann', 'valor': 0.0},
    'bottom': {'type': 'neumann', 'valor': 0.0},
    'top':    {'type': 'neumann', 'valor': 0.0}
}

solver = HeatEquation2D(
    N=150, M=200, T=0.5,
    condiciones_frontera=condiciones
)

soluciones, tiempos = solver.solve(DistribucionGaussiana)
```

### Verificación de Simetría

```python
# La solución debe ser simétrica respecto a x=0.5 e y=0.5
sol_final = soluciones[-1]
N = sol_final.shape[0]

# Simetría en x
diff_x = np.abs(sol_final - np.fliplr(sol_final))
print(f"Error de simetría en x: {diff_x.max():.8f}")

# Simetría en y
diff_y = np.abs(sol_final - np.flipud(sol_final))
print(f"Error de simetría en y: {diff_y.max():.8f}")
```

**Resultado esperado**: Errores \(< 10^{-10}\) (errores de redondeo numérico)

---

## Resumen de Casos de Prueba

| Caso | Condiciones Iniciales | Fronteras | Comportamiento Clave |
|------|--------------------|-----------|----------------------|
| 1. Enfriamiento | Gaussiana | Dirichlet (T=0) | Temperatura → 0 |
| 2. Redistribución | Cuatro Esquinas | Neumann (aislado) | Energía conservada |
| 3. Flujo Constante | Cero | Neumann (flujo) + Dirichlet | Estado estacionario |
| 4. Mixtas | Anillo | Mixtas asimétricas | Flujo direccional |
| 5. Periódico | Gaussiana | Dirichlet oscilante | Onda térmica |
| 6. Convergencia | Gaussiana | Dirichlet | Verificación numérica |
| 7. Simetría | Gaussiana | Neumann simétrico | Simetría preservada |

## Script de Prueba Automatizado

```python
def ejecutar_casos_prueba():
    """Ejecuta todos los casos de prueba y genera reporte"""
    
    casos = [
        ("Enfriamiento", caso1_enfriamiento),
        ("Redistribución", caso2_redistribucion),
        ("Flujo constante", caso3_flujo_constante),
        ("Mixtas", caso4_mixtas),
        ("Periódico", caso5_periodico),
        ("Convergencia", caso6_convergencia),
        ("Simetría", caso7_simetria)
    ]
    
    resultados = []
    
    for nombre, funcion_caso in casos:
        print(f"\n{'='*60}")
        print(f"Ejecutando: {nombre}")
        print('='*60)
        
        try:
            resultado = funcion_caso()
            resultados.append((nombre, "✓ Exitoso", resultado))
            print(f"✓ {nombre} completado exitosamente")
        except Exception as e:
            resultados.append((nombre, "✗ Falló", str(e)))
            print(f"✗ {nombre} falló: {e}")
    
    # Generar reporte
    print("\n" + "="*60)
    print("RESUMEN DE PRUEBAS")
    print("="*60)
    for nombre, estado, info in resultados:
        print(f"{estado} {nombre}")
    
    return resultados

# Ejecutar todas las pruebas
resultados = ejecutar_casos_prueba()
```

---

**Próximo paso**: Revisa [Análisis de Resultados](analisis.md) para técnicas de post-procesamiento y visualización avanzada.