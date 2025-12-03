// Versión paralela (OpenMP) con selección de condición de frontera y condición inicial
// Mantiene esquema Red-Black Gauss-Seidel con tiling (bloques) y validación mínima de entradas

#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <functional>
#include <omp.h>
#include <algorithm> // std::min
#include <string>

using namespace std;

// Clase EcuacionesCalor2D 
// Resuelve la ecuación de calor 2D con condiciones de frontera
// seleccionables y usando Red-Black Gauss-Seidel con tiling + OpenMP
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
    double obtenerValorFrontera(const vector<double>& current, int idx, int vecino, string direccion) const {
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

        // Inicializar coordenadas en paralelo (operación segura y simple)
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < N; i++) {
            xCoord[i] = hx * (i + 1); // evitando bordes
            yCoord[i] = hy * (i + 1); // evitando bordes
        }
    }

    // solveRedBlackStencil (paralela) 
    // Implementa Red-Black Gauss-Seidel con tiling (bloques) y OpenMP.
    // Entrada: uPrev (vector 1D con solución en tiempo n)
    // Salida: vector con u^{n+1} aproximado tras convergencia iterativa
    vector<double> solveRedBlackStencil(const vector<double>& uPrev, int maxIter = 10000, double tol = 1e-10){
        // current: aproximación iterativa, inicializada con uPrev
        vector<double> current = uPrev;
        // inverso de la diagonal (idéntico para todos los nodos internos)
        const double invDiag = 1.0 / (1.0 + 2.0 * (lambdaX + lambdaY));
        // Tamaño de bloque (tiling). Ajustable
        const int BX = 32;
        const int BY = 32;
        // Número de bloques en cada dirección (redondeo hacia arriba)
        const int nBlocksX = (N + BX - 1) / BX;
        const int nBlocksY = (N + BY - 1) / BY;
        // Iteración del solver
        for (int iter = 0; iter < maxIter; iter++) {
            double maxDiff = 0.0;

            // Procesar ambas fases (roja y negra)
            for (int fase = 0; fase < 2; fase++) {
                // Fase roja (i+j par): fase == 0
                // Fase negra (i+j impar): fase == 1
                // Paralelizamos por bloques. Cada hilo procesa ciertos bloques.
                #pragma omp parallel for collapse(2) reduction(max:maxDiff) schedule(static)
                for (int bjy = 0; bjy < nBlocksY; bjy++) {
                    for (int bxi = 0; bxi < nBlocksX; bxi++) {
                        int jStart = bjy * BY;
                        int jEnd   = std::min(N, jStart + BY);
                        int iStart = bxi * BX;
                        int iEnd   = std::min(N, iStart + BX);

                        // Procesar nodos dentro del bloque
                        for (int j = jStart; j < jEnd; j++) {
                            for (int i = iStart; i < iEnd; i++) {
                                // Solo nodos de la fase correspondiente
                                if ((i + j) % 2 != fase) continue;
                                int idx = j * N + i;
                                // Sumas de vecinos en x e y
                                double sumX = 0.0, sumY = 0.0;
                                // vecino izquierda (i-1) o condición de frontera
                                if (i > 0) {
                                    sumX += current[idx - 1];
                                } else {
                                    // primera columna interior -> usar frontera
                                    sumX += obtenerValorFrontera(current, idx, idx + 1, "izq");
                                }
                                // vecino derecha (i+1) o condición de frontera
                                if (i < N - 1) {
                                    sumX += current[idx + 1];
                                } else {
                                    sumX += obtenerValorFrontera(current, idx, idx - 1, "der");
                                }
                                // vecino abajo (j-1) o condición de frontera
                                if (j > 0) {
                                    sumY += current[idx - N];
                                } else {
                                    sumY += obtenerValorFrontera(current, idx, idx + N, "inf");
                                }
                                // vecino arriba (j+1) o condición de frontera
                                if (j < N - 1) {
                                    sumY += current[idx + N];
                                } else {
                                    sumY += obtenerValorFrontera(current, idx, idx - N, "sup");
                                }
                                // Actualización local (solución del stencil)
                                double newValue = (uPrev[idx] + lambdaX * sumX + lambdaY * sumY) * invDiag;
                                // Calcular diferencia para convergencia
                                double diff = fabs(newValue - current[idx]);
                                if (diff > maxDiff) maxDiff = diff;
                                // Guardar nuevo valor
                                current[idx] = newValue;
                            }
                        }
                    }
                }// end parallel (fase roja/negra)
            }
            // Si la máxima diferencia es menor que la tolerancia -> convergencia
            if (maxDiff < tol) break;
        }
        // Devolver aproximación para u^{n+1}
        return current;
    }

    // Método solve
    // Construye la condición inicial evaluando u0Func sobre la malla interior,
    // Recorre los pasos de tiempo, llama al solver iterativo por paso,
    // y guarda snapshots cada saveEvery pasos (incluye t=0 y t=T)
    pair<vector<vector<vector<double>>>, vector<double>> solve(function<double(double,double)> u0Func, int saveEvery){
        // vector 1D que guarda la condición inicial en orden row-major (j*N + i)
        vector<double> initialU(N * N);

        // Rellenar initialU en paralelo
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                initialU[j * N + i] = u0Func(xCoord[i], yCoord[j]);
            }
        }
        // 'soluciones' almacenará matrices N x N para cada snapshot guardado
        vector<vector<vector<double>>> soluciones;
        vector<double> tiempos;

        // Guardar snapshot t=0
        vector<vector<double>> mat0(N, vector<double>(N));
        #pragma omp parallel for collapse(2) schedule(static)
        for (int j = 0; j < N; j++) {
            for (int i = 0; i < N; i++) {
                mat0[j][i] = initialU[j * N + i];
            }
        }
        soluciones.push_back(mat0);
        tiempos.push_back(0.0);

        // uPrev contiene la solución en el paso actual (comienza en t=0)
        vector<double> uPrev = initialU;

        // Bucle temporal: para n = 1..M
        for (int n = 1; n <= M; n++) {
            double tN = n * k;

            // resolver iterativamente el sistema para obtener uNext
            vector<double> uNext = solveRedBlackStencil(uPrev, 1000, 1e-10);

            // guardar snapshot si corresponde
            if (n % saveEvery == 0 || n == M) {
                vector<vector<double>> mat(N, vector<double>(N));
                #pragma omp parallel for collapse(2) schedule(static)
                for (int j = 0; j < N; j++) {
                    for (int i = 0; i < N; i++) {
                        mat[j][i] = uNext[j * N + i];
                    }
                }
                soluciones.push_back(mat);
                tiempos.push_back(tN);
            }
            // Avanzar en el tiempo (mover uNext a uPrev)
            uPrev = std::move(uNext);
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
// Imprime Tmax, Tmin y rango para cada snapshot guardado.
void analizarMinMax(const vector<vector<vector<double>>>& sol, const vector<double>& tiempos)
{
    cout << "\nTiempo      Tmax         Tmin         Rango\n";
    for (size_t k = 0; k < sol.size(); k++) {
        double tmax = -1e300;
        double tmin =  1e300;
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
double pedirNumero(const string &msg){
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
int main(){
    cout << " Simulación Ecuación de Calor 2D (paralela, Red-Black + tiling) \n\n";

    // Parámetros físicos
    double a = 1.0, b = 1.0, c = 1.0; // largo en X y Y, coeficiente (difusividad/velocidad según la formulación)

    // Parámetros numéricos
    int N = 250, M = 8000; // nodos interiores por eje, pasos de tiempo
    double T = 0.4; // tiempo final
    int saveEvery = 2000; // cada cuantos pasos guardar (>=M -> solo t=0 y t=T)

    cout << "Parámetros físicos y numéricos:\n";
    cout << "  a = " << a << ", b = " << b << ", c = " << c << "\n";
    cout << "  N = " << N << ", M = " << M << ", T = " << T << ", saveEvery = " << saveEvery << "\n\n";

    // Lista de condiciones iniciales 
    vector<pair<string,function<double(double,double)>>> testsU0;
    testsU0.push_back(make_pair(string("Distribución Gaussiana"), distribucionGaussiana));
    testsU0.push_back(make_pair(string("Anillo"), anillo));
    testsU0.push_back(make_pair(string("Diagonal"), diagonal));
    testsU0.push_back(make_pair(string("Forma de U"), formaU));
    testsU0.push_back(make_pair(string("Cuatro Esquinas"), cuatroEsquinas));

    // Menu: escoger condición de frontera
    int opcionF;
    double topB = 0.0, botB = 0.0, leftB = 0.0, rightB = 0.0;
    string tipoF;
    cout << "Seleccione condición de frontera:\n";
    cout << " 1. Dirichlet con u=0 en los bordes\n";
    cout << " 2. Dirichlet con temperatura fija\n";
    cout << " 3. Neumann (aislada)\n";
    cout << "Opción: ";
    while (!(cin >> opcionF) || opcionF < 1 || opcionF > 3) {
        cin.clear();
        cin.ignore(10000, '\n');
        cout << "Opción invalida. Intente nuevamente (1-3): ";
    }
    if (opcionF == 1) {
        tipoF = "Dirichlet0";
    }
    else if (opcionF == 2) {
        tipoF = "DirichletCte";
        // Pedir temperaturas de los cuatro bordes (con validación)
        topB   = pedirNumero("Temperatura borde superior: ");
        botB   = pedirNumero("Temperatura borde inferior: ");
        leftB  = pedirNumero("Temperatura borde izquierdo: ");
        rightB = pedirNumero("Temperatura borde derecho: ");
    }
    else { // opcionF == 3
        tipoF = "Neumann";
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

    // Imprimir número de hilos OpenMP que se usarán (sin crear equipo de hilos innecesario)
    int num_procs = 1;
    #pragma omp parallel
    {
        num_procs = omp_get_num_threads();
    }
    cout << "\nOpenMP threads = " << num_procs << "\n";

    // Ejecutar la simulación 
    cout << "\nTEST: " << testsU0[opcionU0 - 1].first << "\n";
    cout << "  Ejecutando con N=" << N << ", M=" << M << ", T=" << T << "\n";

    // Construir el solver con fronteras y parámetros
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
