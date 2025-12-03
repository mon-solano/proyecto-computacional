#include <iostream> 
#include <vector>  
#include <cmath>
#include <chrono>
#include <iomanip>
#include <functional>

using namespace std;

//  Clase EcuacionesCalor2D 
//  Para resolver la ecuación de calor en 2D con condiciones de fronteras seleccionables por el usuario 
class EcuacionesCalor2D {
private:
    // Parámetros de discretización y físicos
    int N, M; // número de nodos interiores por eje (malla N x N), número de pasos de tiempo
    double T, a, b, c; // tiempo final, largo en X y Y del dominio, parámetro físico (velocidad/difusividad según la formulación)

    // Pasos y coeficientes derivados
    double hx, hy, k; // espaciado en X y Y, paso temporal
    double lambdaX; // coeficiente lambda en x
    double lambdaY; // coeficiente lambda en y
    vector<double> xCoord, yCoord; // coordenadas X y Y de los nodos interiores

    // Condiciones de frontera (configurables)
    string tipoFrontera; // "Dirichlet0", "DirichletCte", o "Neumann"
    double bordeSuperior, bordeInferior, bordeIzquierdo, bordeDerecho;

    // Método auxiliar para obtener valor de frontera
    double obtenerValorFrontera(const vector<double>& current, int idx, int vecino, string direccion) {
        if (tipoFrontera == "Dirichlet0") return 0.0;
        if (tipoFrontera == "DirichletCte") {
            if (direccion == "izq") return bordeIzquierdo;
            if (direccion == "der") return bordeDerecho;
            if (direccion == "inf") return bordeInferior;
            if (direccion == "sup") return bordeSuperior;
        }
        // Neumann: usa el vecino opuesto
        return current[vecino];
    }

public:
    // Constructor: recibe N, M, T y opcionales a,b,c y parámetros de frontera
    // Se guardan los parámetros físicos y de entrada y las condiciones de frontera:
    // cadena que indica el tipo elegido, valor top (y = b), valor bottom (y = 0), valor left (x = 0), valor right (x = a)
    EcuacionesCalor2D(int Nparam, int Mparam, double Tparam,
                      double aParam = 1.0, double bParam = 1.0, double cParam = 1.0, string tipoF = "Dirichlet0", 
                      double top = 0.0, double bottom = 0.0, double left = 0.0, double right = 0.0) 
                      : N(Nparam), M(Mparam), T(Tparam), a(aParam), b(bParam), c(cParam),
                        tipoFrontera(tipoF), bordeSuperior(top), bordeInferior(bottom), 
                        bordeIzquierdo(left), bordeDerecho(right){
        // Calcular discretizaciones espaciales
        // Se divide a en N+1 para contar bordes (los nodos interiores son N)
        hx = a / (N + 1);
        hy = b / (N + 1);

        // Paso temporal (k) a partir del tiempo final y número de pasos
        k = T / M;

        // Coeficientes lambda (aparecen tras discretizar las derivadas)
        lambdaX = (c * c) * k / (hx * hx);
        lambdaY = (c * c) * k / (hy * hy);

        // Rellenar vectores de coordenadas (nodos interiores; evitamos fronteras)
        xCoord.resize(N);
        yCoord.resize(N);
        for (int i = 0; i < N; i++) {
        xCoord[i] = hx * (i + 1); // i+1 para evitar el borde x=0
        yCoord[i] = hy * (i + 1); // j+1 para evitar el borde y=0
        }
    }

    //  Método Red-Black Gauss–Seidel (Sin paralelización)
    //  Implementa Red-Black Gauss--Seidel para resolver el sistema lineal que surge del paso implícito de Euler.
    // Entrada: uPrev (vector 1D con solución en tiempo n)
    // Salida: vector con u^{n+1} aproximado tras convergencia del iterativo
    vector<double> solveRedBlackStencil(const vector<double>& uPrev, int maxIter = 10000, double tol = 1e-10) {
        // 'current' contiene la aproximación iterativa; inicializamos con uPrev
        vector<double> current = uPrev;
        // inverso del coeficiente diagonal
        double invDiag = 1.0 / (1.0 + 2.0 * ( lambdaX + lambdaY ));
        // iteración del método hasta maxIter o hasta que maxDiff < tol
        for (int iter = 0; iter < maxIter; iter++) {
            double maxDiff = 0.0;

            // Procesar ambas fases (roja y negra)
            for (int fase = 0; fase < 2; fase++) {
                // Fase roja: (i+j) par (fase==0)
                // Fase negra: (i+j) impar (fase==1)
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < N; i++) {
                        if ((i + j) % 2 == fase) { // condición tablero ajedrez
                            int idx = j * N + i; // índice lineal en vector 1D
                            double sumX = 0.0, sumY = 0.0;

                            // vecino izquierda (i-1)
                            if (i > 0) {
                                sumX += current[idx - 1]; // si existe vecina interior, se usa
                            } else {
                                // Si estamos en la primera columna interior, usamos la frontera según tipo escogido
                                sumX += obtenerValorFrontera(current, idx, idx + 1, "izq");
                            }
                            // vecino derecha (i+1)
                            if (i < N - 1) {
                                sumX += current[idx + 1];
                            } else {
                                // ultima columna interior (i == N-1)
                                sumX += obtenerValorFrontera(current, idx, idx - 1, "der");
                            }
                            // vecino abajo (j-1)
                            if (j > 0) {
                                sumY += current[idx - N];
                            } else {
                                // primera fila interior (j == 0)
                                sumY += obtenerValorFrontera(current, idx, idx + N, "inf");
                            }
                            // vecino arriba (j+1)
                            if (j < N - 1) {
                                sumY += current[idx + N];
                            } else {
                                // ultima fila interior (j == N-1)
                                sumY += obtenerValorFrontera(current, idx, idx - N, "sup");
                            }
                            // Fórmula de actualización (solución del stencil local)
                            double newValue = (uPrev[idx] + lambdaX * sumX + lambdaY * sumY) * invDiag;
                            // Diferencia para criterio de convergencia
                            double diff = fabs(newValue - current[idx]);
                            if (diff > maxDiff) maxDiff = diff;
                            // actualizar valor
                            current[idx] = newValue;
                        }
                    }
                }
            }
            // Si la mayor diferencia es menor que la tolerancia, consideramos convergente
            if (maxDiff < tol) break;
        }
        // Devolvemos la aproximación convergente para u^{n+1}
        return current;
    }

    // Método solve
    // Construye la condición inicial evaluando u0Func sobre la malla interior,
    // Recorre los pasos de tiempo, llama al solver iterativo por paso,
    // y guarda snapshots cada saveEvery pasos (incluye t=0 y t=T)
    pair<vector<vector<vector<double>>>, vector<double>> solve(function<double(double,double)> u0Func, int saveEvery) {
        // vector 1D que guarda la condición inicial en orden row-major (j*N + i)
        vector<double> initialU(N*N);
        for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++)
                initialU[j*N + i] = u0Func(xCoord[i], yCoord[j]);

        // 'soluciones' almacenará matrices N x N para cada snapshot guardado
        vector<vector<vector<double>>> soluciones;
        vector<double> tiempos;

        // Guardar t=0
        vector<vector<double>> mat0(N, vector<double>(N));
        for (int j = 0; j < N; j++)
            for (int i = 0; i < N; i++)
                mat0[j][i] = initialU[j*N + i];
        soluciones.push_back(mat0);
        tiempos.push_back(0.0);

        // uPrev contiene la solución en el paso actual (comienza en t=0)
        vector<double> uPrev = initialU;

        // Bucle temporal: para n = 1..M
        for (int n = 1; n <= M; n++) {
            double tN = n * k; // tiempo actual

            // resolver iterativamente el sistema para obtener uNext
            vector<double> uNext = solveRedBlackStencil(uPrev, 1000, 1e-10);

            // guardar snapshot si corresponde
            if (n % saveEvery == 0 || n == M) {
                vector<vector<double>> mat(N, vector<double>(N));
                for (int j = 0; j < N; j++)
                    for (int i = 0; i < N; i++)
                        mat[j][i] = uNext[j*N + i];
                soluciones.push_back(mat);
                tiempos.push_back(tN);
            }
            // Avanzar en el tiempo (mover uNext a uPrev)
            uPrev = uNext;
        }
        return make_pair(soluciones, tiempos);
    }
}; // fin clase EcuacionesCalor2D

// Condiciones iniciales (funciones que devuelven u0(x,y))
// Gaussiana centrada en (0.5,0.5) con amplitud 20
double distribucionGaussiana(double x, double y) {
    return 20.0 * exp(-50.0*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
}
// Anillo centrado en (0.5,0.5) con radio ~0.3
double anillo(double x, double y) {
    double r = sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
    return 20.0 * exp(-200.0*(r-0.3)*(r-0.3));
}
// Distribución sobre la diagonal x=y
double diagonal(double x, double y) {
    return 20.0 * exp(-100.0*(x-y)*(x-y));
}
// Forma en U (tres picos)
double formaU(double x, double y) {
    double l = exp(-200.0*((x-0.25)*(x-0.25)+(y-0.5)*(y-0.5)));
    double r = exp(-200.0*((x-0.75)*(x-0.75)+(y-0.5)*(y-0.5)));
    double b = exp(-200.0*((x-0.5)*(x-0.5)+(y-0.25)*(y-0.25)));
    return 20.0*(l + r + b);
}
// Cuatro picos en las esquinas interiores
double cuatroEsquinas(double x, double y) {
    double tl = exp(-200.0*((x-0.2)*(x-0.2)+(y-0.8)*(y-0.8)));
    double tr = exp(-200.0*((x-0.8)*(x-0.8)+(y-0.8)*(y-0.8)));
    double bl = exp(-200.0*((x-0.2)*(x-0.2)+(y-0.2)*(y-0.2)));
    double br = exp(-200.0*((x-0.8)*(x-0.8)+(y-0.2)*(y-0.2)));
    return 20.0*(tl + tr + bl + br);
}

// Función auxiliar: analizarMinMax
// Imprime Tmax, Tmin y rango para cada snapshot guardado
void analizarMinMax(const vector<vector<vector<double>>>& sol, const vector<double>& tiempos)
{
    cout << "\nTiempo      Tmax         Tmin         Rango\n";
    for (size_t k = 0; k < sol.size(); k++) {
        double tmax = -1e300, tmin = 1e300;
        // recorrer cada fila (cada matriz) y encontrar max/min
        for (const auto& fila : sol[k]) {
            for (double v : fila) {
                if (v > tmax) tmax = v;
                if (v < tmin) tmin = v;
            }
        }
        cout << fixed << setprecision(4);
        cout << tiempos[k] << "      " 
             << tmax << "       "
             << tmin << "       "
             << (tmax - tmin) << "\n";
    }
}

// Funcion para pedir número
double pedirNumero(const string &msg) {
    double val;
    while (true) {
        cout << msg;
        cin >> val;
        if (!cin.fail()) break;      // La entrada es válida (numérica)
        cin.clear();                 // Limpia el estado de error
        cin.ignore(10000, '\n');     // Descarta basura en el buffer
        cout << "Entrada invalida. Ingrese un numero.\n";
    }
    return val;
}

// main: orquesta la simulación, permite seleccionar frontera
int main() {
    cout << " Simulación Ecuación de Calor 2D \n";

    // Parámetros físicos
    double a = 1.0, b = 1.0, c = 1.0; // largo en X y Y, coeficiente (difusividad/velocidad según la formulación)

    // Parámetros numéricos
    int N = 50, M = 500; // nodos interiores por eje, pasos de tiempo
    double T = 0.4; // tiempo final
    int saveEvery = 2000; // cada cuantos pasos guardar (>=M -> solo t=0 y t=T)

    cout << "Parámetros físicos y numéricos:\n";
    cout << "  a = " << a << ", b = " << b << ", c = " << c << "\n";
    cout << "  N = " << N << ", M = " << M << ", T = " << T << ", saveEvery = " << saveEvery << "\n";

    // Lista de condiciones iniciales
    vector<pair<string,function<double(double,double)>>> testsU0 = {
        {"Distribución Gaussiana", distribucionGaussiana},
        {"Anillo", anillo},
        {"Diagonal", diagonal},
        {"Forma de U", formaU},
        {"Cuatro Esquinas", cuatroEsquinas}
    };

    // Menu: escoger condición de frontera
    int opcion;
    double topB=0, botB=0, leftB=0, rightB=0;
    string tipoF;
    cout << "\nSeleccione condición de frontera:\n";
    cout << "1. Dirichlet con u=0 en los bordes\n";
    cout << "2. Dirichlet con temperatura fija\n";
    cout << "3. Neumann (aislada)\n";
    cout << "Opción: ";
    while (!(cin >> opcion) || opcion < 1 || opcion > 3) {
        cin.clear();
        cin.ignore(10000, '\n');
        cout << "Opción invalida. Intente nuevamente (1-3): ";
    }
    if(opcion == 1) {
        tipoF = "Dirichlet0"; // frontera homogénea
    }
    else if(opcion == 2) {
        tipoF = "DirichletCte"; // frontera con valores fijos
        topB   = pedirNumero("Temperatura borde superior: ");
        botB   = pedirNumero("Temperatura borde inferior: ");
        leftB  = pedirNumero("Temperatura borde izquierdo: ");
        rightB = pedirNumero("Temperatura borde derecho: ");
    }
    else if(opcion == 3) {
        tipoF = "Neumann"; // frontera aislada
    }

    // Menú para condición inicial
    int opcionU0;
    cout << "\nSeleccione condición inicial:\n";
    for (int i = 0; i < (int)testsU0.size(); i++) {
        cout << " " << (i + 1) << ". " << testsU0[i].first << "\n";
    }
    cout << "Opción: ";
    while (!(cin >> opcionU0) || opcionU0 < 1 || opcionU0 > (int)testsU0.size()) {
        cin.clear();
        cin.ignore(10000, '\n');
        cout << "Opción invalida. Intente nuevamente: ";
    }

    // Seleccionar la función inicial elegida
    function<double(double,double)> u0Func = testsU0[opcionU0 - 1].second;

    // Ejecutar el test eligido
    cout << "\nTEST: " << testsU0[opcionU0 - 1].first << "\n";
    cout << "  Ejecutando con N=" << N << ", M=" << M << ", T=" << T << "\n";

    // Construir el solver con parámetros y condiciones de frontera        
    EcuacionesCalor2D solver(N, M, T, a, b, c, tipoF, topB, botB, leftB, rightB);

    // Tiempo del test
    auto start = chrono::high_resolution_clock::now();

    // Ejecutar la simulación
    pair<vector<vector<vector<double>>>, vector<double>> solTiempo = solver.solve(u0Func, saveEvery);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    cout << "  Tiempo: " << fixed << setprecision(2) << elapsed.count() << " s\n";

    // Analizar resultados guardados 
    analizarMinMax(solTiempo.first, solTiempo.second);

    cout << "\n Simulación completada.\n";
    return 0;
}
