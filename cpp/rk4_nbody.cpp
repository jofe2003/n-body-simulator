/*
 * Simulador gravitacional de N corpos integrado com Runge-Kutta classico de 4a ordem.
 *
 * Uso:  ./rk4_nbody [entrada] [saida]        (padrao: data.txt e output.txt)
 *
 * Formato da entrada (linhas iniciadas com '#' sao comentarios):
 *     n_corpos  dt  n_passos  salvar_a_cada
 *     x y z vx vy vz massa          <- uma linha por corpo (unidades SI)
 *
 * Saida: uma linha por instante salvo com "x y z" de cada corpo, em sequencia.
 */
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

const double G = 6.67430e-11;

/* Estado do sistema: para cada corpo, (x, y, z, vx, vy, vz) em sequencia */
typedef vector<double> State;

/* Calcula dY/dt: derivada da posicao = velocidade; derivada da velocidade = aceleracao */
void derivative(const State &Y, const vector<double> &mass, State &dY) {
    const size_t n = mass.size();

    for (size_t i = 0; i < n; i++) {
        for (int k = 0; k < 3; k++) {
            dY[6*i + k]     = Y[6*i + 3 + k];
            dY[6*i + 3 + k] = 0.0;
        }
    }

    /* Cada par (i, j) e calculado uma unica vez (3a lei de Newton) */
    for (size_t i = 0; i < n; i++) {
        for (size_t j = i + 1; j < n; j++) {
            double r[3], r2 = 0.0;
            for (int k = 0; k < 3; k++) {
                r[k] = Y[6*j + k] - Y[6*i + k];
                r2 += r[k] * r[k];
            }
            const double invR3 = 1.0 / (r2 * sqrt(r2));

            for (int k = 0; k < 3; k++) {
                dY[6*i + 3 + k] += G * mass[j] * r[k] * invR3;  /* i e puxado para j */
                dY[6*j + 3 + k] -= G * mass[i] * r[k] * invR3;  /* j e puxado para i */
            }
        }
    }
}

/* Um passo de Runge-Kutta 4: Y(t+dt) = Y(t) + dt*(k1 + 2k2 + 2k3 + k4)/6 */
void rungeKutta4(State &Y, const vector<double> &mass, double dt) {
    const size_t m = Y.size();
    State k1(m), k2(m), k3(m), k4(m), tmp(m);

    derivative(Y, mass, k1);
    for (size_t q = 0; q < m; q++) tmp[q] = Y[q] + 0.5 * dt * k1[q];
    derivative(tmp, mass, k2);
    for (size_t q = 0; q < m; q++) tmp[q] = Y[q] + 0.5 * dt * k2[q];
    derivative(tmp, mass, k3);
    for (size_t q = 0; q < m; q++) tmp[q] = Y[q] + dt * k3[q];
    derivative(tmp, mass, k4);

    for (size_t q = 0; q < m; q++)
        Y[q] += dt * (k1[q] + 2.0*k2[q] + 2.0*k3[q] + k4[q]) / 6.0;
}

/* Le o arquivo de entrada ignorando comentarios */
bool readInput(const string &path, State &Y, vector<double> &mass,
               double &dt, long &steps, long &saveEvery) {
    ifstream file(path);
    if (!file) return false;

    stringstream clean;
    string line;
    while (getline(file, line)) {
        size_t hash = line.find('#');
        if (hash != string::npos) line.erase(hash);
        clean << line << '\n';
    }

    int n;
    if (!(clean >> n >> dt >> steps >> saveEvery) || n < 1 || saveEvery < 1) return false;

    Y.assign(6 * n, 0.0);
    mass.assign(n, 0.0);
    for (int i = 0; i < n; i++) {
        for (int k = 0; k < 6; k++)
            if (!(clean >> Y[6*i + k])) return false;
        if (!(clean >> mass[i])) return false;
    }
    return true;
}

void writePositions(ofstream &out, const State &Y, size_t n) {
    for (size_t i = 0; i < n; i++)
        for (int k = 0; k < 3; k++)
            out << Y[6*i + k] << ' ';
    out << '\n';
}

int main(int argc, char **argv) {
    const string inputPath  = argc > 1 ? argv[1] : "data.txt";
    const string outputPath = argc > 2 ? argv[2] : "output.txt";

    State Y;
    vector<double> mass;
    double dt;
    long steps, saveEvery;

    if (!readInput(inputPath, Y, mass, dt, steps, saveEvery)) {
        cerr << "Erro ao ler a entrada '" << inputPath << "'\n";
        return 1;
    }

    ofstream output(outputPath);
    if (!output) {
        cerr << "Erro ao criar a saida '" << outputPath << "'\n";
        return 1;
    }
    output << scientific << setprecision(15);

    auto start = chrono::high_resolution_clock::now();

    writePositions(output, Y, mass.size());
    for (long i = 1; i <= steps; i++) {
        rungeKutta4(Y, mass, dt);
        if (i % saveEvery == 0) writePositions(output, Y, mass.size());
    }

    chrono::duration<double> elapsed = chrono::high_resolution_clock::now() - start;
    cout << mass.size() << " corpos, " << steps << " passos em "
         << elapsed.count() << " s -> " << outputPath << endl;
    return 0;
}
