/*
 * Simulador gravitacional de N corpos integrado com o metodo de Euler semi-implicito.
 *
 * Uso:  ./euler_nbody [entrada] [saida]       (padrao: data.txt e output.txt)
 * Formato da entrada: o mesmo de rk4_nbody.cpp.
 */
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <sstream>
#include <string>
using namespace std;

struct bodyInfo {
    double position[3] {0};
    double velocity[3] {0};
    double mass {0};
    double *distance {0};
    double force[3] {0};
    double acceleration[3] {0};
};

void calculateDistance(bodyInfo *body, int numberOfBodies);
void calculateForce(bodyInfo *body, int numberOfBodies);
void calculateAcceleration(bodyInfo *body, int numberOfBodies);
void update(bodyInfo *body, double timeInterval, int numberOfBodies);

int main(int argc, char **argv) {
    /* Uso: ./euler_nbody [entrada] [saida]   (padrao: data.txt e output.txt) */
    const string inputPath  = argc > 1 ? argv[1] : "data.txt";
    const string outputPath = argc > 2 ? argv[2] : "output.txt";

    int numberOfBodies {0};
    long steps {0}, saveEvery {1};
    double number {0}, timeInterval {0};

    /* Remove comentarios (#) do arquivo de entrada */
    ifstream file (inputPath);
    if (!file) {
        cerr << "Erro ao ler a entrada '" << inputPath << "'" << endl;
        return 1;
    }
    stringstream data;
    string line;
    while (getline(file, line)) {
        size_t hash = line.find('#');
        if (hash != string::npos) line.erase(hash);
        data << line << '\n';
    }

    /* Os primeiros dados a serem lidos sao o numero de corpos, o intervalo de
       tempo (delta T), o numero de passos e a cada quantos passos salvar */
    data >> numberOfBodies >> timeInterval >> steps >> saveEvery;
    if (numberOfBodies < 1 || saveEvery < 1) {
        cerr << "Cabecalho invalido em '" << inputPath << "'" << endl;
        return 1;
    }

    /* Creates an array of bodyInfo of the size of numberOfBodies */
    bodyInfo* body = new bodyInfo[numberOfBodies];

    /* Creates and initializes de distance double for each body */
    for (int i = 0; i < numberOfBodies; ++i) {
        body[i].distance = new double[numberOfBodies];

        for (int j = 0; j < numberOfBodies; ++j) {
            body[i].distance[j] = 0;
        }
    }

    /* Stores the initial data for position, velocity and mass of each body */
    for(int i = 0; i < numberOfBodies; i++) {
        for(int j = 0; j < 3; j++) {
            data >> number;
            body[i].position[j] = number;
        }

        for(int j = 0; j < 3; j++) {
            data >> number;
            body[i].velocity[j] = number;
        }

        data >> number;
        body[i].mass = number;
    }

    ofstream output (outputPath);
    output << scientific << setprecision(15);

    auto start_time = std::chrono::high_resolution_clock::now();

    for (long i = 0; i <= steps; i++) {
        if (i % saveEvery == 0) {
            for (int b = 0; b < numberOfBodies; b++) {
                for (int k = 0; k < 3; k++) {
                    output << body[b].position[k] << " ";
                }
            }
            output << "\n";
        }

        if (i == steps) break;

        calculateDistance(body, numberOfBodies);

        calculateForce(body, numberOfBodies);

        calculateAcceleration(body, numberOfBodies);

        update(body, timeInterval, numberOfBodies);
    }

    for (int i = 0; i < numberOfBodies; ++i) {
        delete[] body[i].distance;
    }
    delete[] body;
    output.close();

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
    std::cout << numberOfBodies << " corpos, " << steps << " passos em "
              << duration.count() << " s -> " << outputPath << std::endl;

    return 0;
}

/*calcula a distancia entre os corpos*/
void calculateDistance(bodyInfo *body, int numberOfBodies) {
    double x, y, z;

    /* Here, one body [i] will be selected to have its distance calculated compared
       to every other body in the system */
    for (int i = 0; i < numberOfBodies; i++) {
        /* This is where body [i] will have its distance calculated in respect
            to all the other [j] bodies */
        for (int j = i + 1; j < numberOfBodies; j++) {
            x = body[i].position[0] - body[j].position[0];
            y = body[i].position[1] - body[j].position[1];
            z = body[i].position[2] - body[j].position[2];
            
            body[i].distance[j] = sqrt((x*x) + (y*y) + (z*z));
            body[j].distance[i] = body[i].distance[j];
        }
    }

    return;       
}

void calculateForce(bodyInfo *body, int numberOfBodies) {
    double mass1, mass2, distanceCubed, distanceSubtracted;

    /* Here, one body [i] will be selected to have its distance calculated compared
       to every other body in the system */
    for (int i = 0; i < numberOfBodies; i++) {
        /* Iteration for each on of the coordinates (x, y, z) */
        for (int k = 0; k < 3; k++) {
            /* This is where body [i] will have its distance calculated in respect
               to all the other [j] bodies */
            body[i].force[k] = 0;
            for (int j = 0; j < numberOfBodies; j++) {
                if (i != j) {
                    mass1 = body[i].mass;
                    mass2 = body[j].mass;
                    distanceCubed = pow(body[i].distance[j], 3);
                    distanceSubtracted = body[i].position[k] - body[j].position[k];

                    body[i].force[k] += -6.6743e-11*((mass1*mass2)/distanceCubed)*distanceSubtracted;
                }
            }
        }
    }    

    return;                          
}

void calculateAcceleration(bodyInfo *body, int numberOfBodies) {
    /* Calculate the acceleration of everybody on x, y and z */
    for (int i = 0; i < numberOfBodies; i++) {
        for (int k = 0; k < 3; k++) {
            body[i].acceleration[k] = body[i].force[k] / body[i].mass;
        }
    }
    return;   
}

void update(bodyInfo *body, double timeInterval, int numberOfBodies) {
    /* Updates the velocity and position of each body using the semi-implicit
       (symplectic) Euler method: the new velocity is used to update the position */
    for (int i = 0; i < numberOfBodies; i++) {
        for (int k = 0; k < 3; k++) {
            body[i].velocity[k] = body[i].velocity[k] + (timeInterval * body[i].acceleration[k]);
            body[i].position[k] = body[i].position[k] + (timeInterval * body[i].velocity[k]);
        }
    }
    return;   
}