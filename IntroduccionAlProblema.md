# **ECUACIÓN DE CALOR EN DOS DIMENSIONES**

La idea del proyecto corresponde a resolver de manera numérica la ecuación de calor en dos dimensiones.  
Suponga que u = u (x, y, t) es una variable escalar que define la temperatura de una región de dos dimensiones en el plano Cartesiano (x, y) como función del tiempo. Bajo condiciones ideales (sin fuentes de energía externas, capacidad calórica uniforme, aislamiento perfecto), la ecuación de movimiento está dada por:  

$$
\frac{\partial u}{\partial t} = c^2 \left[ \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right],
$$

Para alguna región acotada $ x \in [0,a] $ y $ y \in [0,b] $. Las condiciones de frontera pueden cambiar y la dinámica está determinada por éstas y por las condiciones iniciales.  

**Objetivos del proyecto:**

- Encontrar una metodología numérica para resolver el sistema de ecuaciones  
- Implementar la solución en **Python**.  
- Implementar la solución en **C++**.  
- Utilice distintas condiciones iniciales.  
- En un gráfico de dos dimensiones con mapa de colores, encuentre una forma de visualizar la solución a través del tiempo. 
- Encontrar una forma de acelerar la aplicación utilizando paralelismo de memoria compartida.  

