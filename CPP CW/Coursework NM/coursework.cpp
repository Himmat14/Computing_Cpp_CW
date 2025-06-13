// C++ Code for Part II of the Coursework
// Numerical Simulation of Vessel Dynamics

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

struct Vessel {
    double mass, area, flow_rate_coeff, y0, v0, h0;
    double y, v, h;

    Vessel(double m, double a, double k, double y_init, double v_init, double h_init)
            : mass(m), area(a), flow_rate_coeff(k), y0(y_init), v0(v_init), h0(h_init), y(y_init), v(v_init), h(h_init) {}
};

void readParameters(const std::string &filename, int &scheme, double &T, double &dt, double &g, double &rho, double &Cd, std::vector<Vessel> &vessels) {
    std::ifstream infile(filename);
    if (!infile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::istringstream iss(line);
        if (scheme == -1) { // Read first line for global parameters
            iss >> scheme >> T >> dt >> g >> rho >> Cd;
        } else { // Read vessel parameters
            double m, a, k, y_init, v_init, h_init;
            iss >> m >> a >> k >> y_init >> v_init >> h_init;
            vessels.emplace_back(m, a, k, y_init, v_init, h_init);
        }
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

void simulateEuler(double T, double dt, double g, double rho, double Cd, std::vector<Vessel> &vessels, std::vector<std::vector<double>> &results) {
    size_t steps = static_cast<size_t>(T / dt);
    for (size_t step = 0; step <= steps; ++step) {
        std::vector<double> snapshot;
        for (auto &v : vessels) {
            double Fw = v.mass * g + rho * g * v.area * v.h;
            double Fb = -rho * g * v.area * (v.h0 + v.y - v.h);
            double Fd = -0.5 * rho * v.area * Cd * std::abs(v.v) * v.v;
            double total_mass = v.mass + rho * v.area * v.h;

            double a = (Fw + Fb + Fd) / total_mass;

            v.v += a * dt;
            v.y += v.v * dt;
            v.h += (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y - v.h) * dt;

            snapshot.push_back(v.y);
            snapshot.push_back(v.v);
            snapshot.push_back(v.h);
        }
        results.push_back(snapshot);
    }
}

void simulateRK4(double T, double dt, double g, double rho, double Cd, std::vector<Vessel> &vessels, std::vector<std::vector<double>> &results) {
    size_t steps = static_cast<size_t>(T / dt);
    for (size_t step = 0; step <= steps; ++step) {
        std::vector<double> snapshot;
        for (auto &v : vessels) {
            auto dynamics = [&](double y, double v, double h, const Vessel &vessel) {
                double Fw = vessel.mass * g + rho * g * vessel.area * h;
                double Fb = -rho * g * vessel.area * (vessel.h0 + y - h);
                double Fd = -0.5 * rho * vessel.area * Cd * std::abs(v) * v;
                double total_mass = vessel.mass + rho * vessel.area * h;
                double a = (Fw + Fb + Fd) / total_mass;
                double dhdt = (vessel.flow_rate_coeff / (rho * vessel.area)) * (vessel.h0 + y - h);
                return std::make_tuple(a, dhdt);
            };

            // RK4 Steps
            double k1_y = v.v;
            double k1_v, k1_h;
            std::tie(k1_v, k1_h) = dynamics(v.y, v.v, v.h, v);

            double k2_y = v.v + dt / 2 * k1_v;
            double k2_v, k2_h;
            std::tie(k2_v, k2_h) = dynamics(v.y + dt / 2 * k1_y, v.v + dt / 2 * k1_v, v.h + dt / 2 * k1_h, v);

            double k3_y = v.v + dt / 2 * k2_v;
            double k3_v, k3_h;
            std::tie(k3_v, k3_h) = dynamics(v.y + dt / 2 * k2_y, v.v + dt / 2 * k2_v, v.h + dt / 2 * k2_h, v);

            double k4_y = v.v + dt * k3_v;
            double k4_v, k4_h;
            std::tie(k4_v, k4_h) = dynamics(v.y + dt * k3_y, v.v + dt * k3_v, v.h + dt * k3_h, v);

            v.y += dt / 6 * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
            v.v += dt / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
            v.h += dt / 6 * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);

            snapshot.push_back(v.y);
            snapshot.push_back(v.v);
            snapshot.push_back(v.h);
        }
        results.push_back(snapshot);
    }
}

int main() {
    int scheme = -1;
    double T, dt, g, rho, Cd;
    std::vector<Vessel> vessels;
    readParameters("test_data.txt", scheme, T, dt, g, rho, Cd, vessels);

    std::vector<std::vector<double>> results;
    if (scheme == 0) {
        simulateEuler(T, dt, g, rho, Cd, vessels, results);
    } else if (scheme == 1) {
        simulateRK4(T, dt, g, rho, Cd, vessels, results);
    } else {
        std::cerr << "Error: Unknown scheme selected." << std::endl;
        return 1;
    }

    writeOutput("output.txt", results, dt);

    return 0;
}
