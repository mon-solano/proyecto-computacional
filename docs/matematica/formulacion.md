# Formulación Matemática

## La Ecuación de Calor

La ecuación de calor en dos dimensiones es una ecuación diferencial parcial parabólica que describe la distribución de temperatura en una región espacial a lo largo del tiempo.

### Ecuación Diferencial

La forma general de la ecuación de calor en 2D es:

\[
\frac{\partial u}{\partial t} = c^2 \left[\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right] 
\]

donde:

- \(u(x, y, t)\): temperatura en el punto \((x, y)\) en el tiempo \(t\)
- \(c^2\): coeficiente de difusión térmica

### Dominio y Condiciones

El problema se define en un dominio rectangular:

\[
D = [0, a] \times [0, b] \subset \mathbb{R}^2
\]

en un intervalo temporal:

\[
t \in [0, T]
\]

## Condiciones Iniciales

La distribución de temperatura en \(t = 0\) debe especificarse:

\[
u(x, y, 0) = u_0(x, y) \quad \forall (x, y) \in D
\]

Esta función \(u_0(x, y)\) determina cómo comienza el proceso de difusión térmica.

## Condiciones de Frontera

Las condiciones de frontera describen cómo interactúa la región con su entorno. Se implementaron dos tipos principales:

### Condiciones de Dirichlet

Describen una **temperatura constante** en la frontera:

\[
u(x, y, t) = g(x, y, t) \quad \text{en } \partial D
\]

**Interpretación física**: La frontera se mantiene a una temperatura conocida y constante.

### Condiciones de Neumann

Especifican el **flujo de calor** en la frontera:

\[
\frac{\partial u}{\partial n} = g(x, y, t) \quad \text{en } \partial D
\]

**Interpretación física**: Se prescribe el flujo térmico a través de la frontera:

- Si g > 0: el calor entra.
- Si g < 0: el calor sale.
- Si g = 0: el sistema es aislado.

---

**Ahora**: En la sección [Metodología Numérica](metodologia.md) se puede ver cómo discretizamos y resolvemos esta ecuación numéricamente.