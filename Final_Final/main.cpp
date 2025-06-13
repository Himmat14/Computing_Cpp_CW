// Main Program File: simulation_main.cpp
#include "E:\Y2\computing CW\Final_Final\cmake-build-debug\Simulation.h"
#include "E:\Y2\computing CW\Final_Final\cmake-build-debug\EulerSimulation.h"
#include "E:\Y2\computing CW\Final_Final\cmake-build-debug\RK4Simulation.h"
#include "E:\Y2\computing CW\Final_Final\cmake-build-debug\Vessel.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

void readParameters(const std::string &filename, int &scheme, double &T, double &dt, double &g, double &rho, double &Cd, std::vector<Vessel> &vessels) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    bool readGlobalParams = false;

    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);

        if (!readGlobalParams) {
            if (!(iss >> scheme >> T >> dt >> g >> rho >> Cd)) {
                std::cerr << "Error: Invalid global parameters format in input file." << std::endl;
                exit(1);
            }
            readGlobalParams = true;
        } else {
            double m, a, k, y_init, v_init, h_init;
            if (!(iss >> m >> a >> k >> y_init >> v_init >> h_init)) {
                std::cerr << "Error: Invalid vessel parameters format in input file." << std::endl;
                exit(1);
            }
            vessels.push_back(Vessel(m, a, k, y_init, v_init, h_init));
        }
    }

    if (!readGlobalParams) {
        std::cerr << "Error: Global parameters not found in input file." << std::endl;
        exit(1);
    }

    if (vessels.empty()) {
        std::cerr << "Error: No vessel data found in input file." << std::endl;
        exit(1);
    }
}

void writeOutput(const std::string &filename, const std::vector<std::vector<double>> &results, double dt) {
    std::ofstream outfile(filename);
    if (!outfile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }

    for (size_t i = 0; i < results.size(); ++i) {
        outfile << std::fixed;
        outfile << i * dt;
        for (const auto &val : results[i]) {
            outfile << " " << val;
        }
        outfile << "\n";
    }
}

int main() {
    int scheme = -1;
    double T, dt, g, rho, Cd;
    std::vector<Vessel> vessels;

    readParameters("test_data.txt", scheme, T, dt, g, rho, Cd, vessels);

    Simulation *simulation;
    if (scheme == 0) {
        simulation = new EulerSimulation();
    } else if (scheme == 1) {
        simulation = new RK4Simulation();
    } else {
        std::cerr << "Error: Unknown scheme selected." << std::endl;
        return 1;
    }

    simulation->setParameters(T, dt, g, rho, Cd);
    simulation->runSimulation(vessels);
    writeOutput("output.txt", simulation->getResults(), dt);

    delete simulation;

    return 0;
}
