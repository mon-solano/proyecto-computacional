# Ecuación de Calor en Dos Dimensiones

## Bienvenido

Este proyecto presenta una solución numérica completa para la **ecuación de calor en dos dimensiones**, implementando métodos de diferencias finitas con el **esquema de Euler Implícito** (Backward Euler).

## Visión General del Proyecto

La ecuación de calor describe cómo se distribuye la temperatura en una región a lo largo del tiempo. En este proyecto, resolvemos:

\[
\frac{\partial u}{\partial t} = c^2 \left[\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2}\right]
\]

donde:

- \(u = u(x, y, t)\): temperatura en el punto \((x, y)\) en el tiempo \(t\)
- \(c\): coeficiente de difusión térmica
- \(x \in [0, a]\), \(y \in [0, b]\): región espacial bidimensional

## Características del proyecto

- Implementación en **Python**
- Animación de la **evolución temporal** para diversas condiciones iniciales y de fronter
- Implementación en **C++**
- Aceleración de la implementación mediante **paralelismo de memoria compartida**


