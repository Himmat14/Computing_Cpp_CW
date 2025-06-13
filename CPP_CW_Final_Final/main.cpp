// C++ Code for Part II of the Coursework
// Numerical Simulation of Vessel Dynamics using OOP

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <sstream>

// Abstract Base Class for Simulation
// Provides an interface for different simulation methods (e.g., Euler, RK4)
class Simulation {
protected:
    double T, dt, g, rho, Cd; // Simulation parameters
    std::vector<std::vector<double>> results; // Store simulation results

public:
    // Pure virtual function to be implemented by derived classes
    virtual void runSimulation(std::vector<class Vessel> &vessels) = 0;

    // Set simulation parameters
    void setParameters(double T, double dt, double g, double rho, double Cd) {
        this->T = T;
        this->dt = dt;
        this->g = g;
        this->rho = rho;
        this->Cd = Cd;
    }

    // Retrieve the simulation results
    const std::vector<std::vector<double>> &getResults() const {
        return results;
    }
};

// Class for Vessel Properties
// Encapsulates the properties and state of a vessel
class Vessel {
public:
    double mass, area, flow_rate_coeff, y0, v0, h0; // Initial properties
    double y, v, h; // Current state (position, velocity, water level)
    double dt; // Time step for updates

    // Constructor to initialize the vessel's properties and state
    Vessel(double m, double a, double k, double y_init, double v_init, double h_init)
            : mass(m), area(a), flow_rate_coeff(k), y0(y_init), v0(v_init), h0(h_init), y(y_init), v(v_init), h(h_init), dt(0) {}

    // Update the vessel's state using acceleration and rate of water level change
    virtual void updateState(double a, double dhdt) {
        v += a * dt; // Update velocity
        y += v * dt; // Update position
        h += dhdt * dt; // Update water level
    }
};

// Derived Class for Explicit Euler Simulation
// Implements the Euler method for numerical integration
class EulerSimulation : public Simulation {
public:
    void runSimulation(std::vector<Vessel> &vessels) override {
        // Assign time step to each vessel
        for (auto &v : vessels) {
            v.dt = dt;
        }

        // Number of time steps
        size_t steps = static_cast<size_t>(T / dt);

        // Simulation loop
        for (size_t step = 0; step <= steps; ++step) {
            std::vector<double> snapshot;
            for (auto &v : vessels) {
                // Calculate forces and acceleration
                double Fw = v.mass * g + rho * g * v.area * v.h; // Weight
                double Fb = -rho * g * v.area * (v.h0 + v.y - v.h); // Buoyancy
                double Fd = -0.5 * rho * v.area * Cd * std::abs(v.v) * v.v; // Drag
                double total_mass = v.mass + rho * v.area * v.h; // Effective mass

                double a = (Fw + Fb + Fd) / total_mass; // Acceleration
                double dhdt = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y - v.h); // Rate of water level change

                // Update vessel state
                v.updateState(a, dhdt);

                // Save current state to snapshot
                snapshot.push_back(v.y);
                snapshot.push_back(v.v);
                snapshot.push_back(v.h);
            }
            results.push_back(snapshot); // Save snapshot to results
        }
    }
};

// Derived Class for Runge-Kutta Simulation
// Implements the 4th-order Runge-Kutta method for numerical integration
class RK4Simulation : public Simulation {
public:
    void runSimulation(std::vector<Vessel> &vessels) override {
        // Assign time step to each vessel
        for (auto &v : vessels) {
            v.dt = dt;
        }

        // Number of time steps
        size_t steps = static_cast<size_t>(T / dt);

        // Simulation loop
        for (size_t step = 0; step <= steps; ++step) {
            std::vector<double> snapshot;
            for (auto &v : vessels) {
                // Lambda function to calculate dynamics (acceleration and dh/dt)
                auto dynamics = [&](double y, double v, double h, const Vessel &vessel) {
                    double Fw = vessel.mass * g + rho * g * vessel.area * h; // Weight
                    double Fb = -rho * g * vessel.area * (vessel.h0 + y - h); // Buoyancy
                    double Fd = -0.5 * rho * vessel.area * Cd * std::abs(v) * v; // Drag
                    double total_mass = vessel.mass + rho * vessel.area * h; // Effective mass
                    double a = (Fw + Fb + Fd) / total_mass; // Acceleration
                    double dhdt = (vessel.flow_rate_coeff / (rho * vessel.area)) * (vessel.h0 + y - h); // Rate of water level change
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

                // Combine RK4 steps to update state
                double a = dt / 6 * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
                double dhdt = dt / 6 * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);

                v.updateState(a, dhdt);

                // Save current state to snapshot
                snapshot.push_back(v.y);
                snapshot.push_back(v.v);
                snapshot.push_back(v.h);
            }
            results.push_back(snapshot); // Save snapshot to results
        }
    }
};

// Helper Functions for File I/O
// Reads simulation parameters and vessel data from a file
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
            vessels.emplace_back(m, a, k, y_init, v_init, h_init);
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

// Writes simulation results to a file
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
    int scheme = -1; // Simulation scheme: 0 = Euler, 1 = RK4
    double T, dt, g, rho, Cd; // Global parameters
    std::vector<Vessel> vessels; // List of vessels

    // Read parameters from file
    readParameters("E:\\Y2\\computing CW\\CPP_CW_Final_Final\\test_data.txt", scheme, T, dt, g, rho, Cd, vessels);

    // Create appropriate simulation object
    Simulation *simulation;
    if (scheme == 0) {
        simulation = new EulerSimulation();
    } else if (scheme == 1) {
        simulation = new RK4Simulation();
    } else {
        std::cerr << "Error: Unknown scheme selected." << std::endl;
        return 1;
    }

    // Set simulation parameters and run
    simulation->setParameters(T, dt, g, rho, Cd);
    simulation->runSimulation(vessels);

    // Write results to file
    writeOutput("output.txt", simulation->getResults(), dt);

    // Clean up
    delete simulation;

    return 0;
}
