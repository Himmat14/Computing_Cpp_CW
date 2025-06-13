// EulerSimulation.h
#ifndef EULERSIMULATION_H
#define EULERSIMULATION_H

#include "Simulation.h"
#include <cmath> // For std::abs

class EulerSimulation : public Simulation {
public:
    void runSimulation(std::vector<Vessel> &vessels) override {
        for (size_t i = 0; i < vessels.size(); ++i) {
            vessels[i].dt = dt;
        }

        size_t steps = static_cast<size_t>(T / dt);

        for (size_t step = 0; step <= steps; ++step) {
            std::vector<double> currentResults;
            for (size_t i = 0; i < vessels.size(); ++i) {
                Vessel &v = vessels[i];
                double Fw = v.mass * g + rho * g * v.area * v.h;
                double Fb = -rho * g * v.area * (v.h0 + v.y - v.h);
                double Fd = -0.5 * rho * v.area * Cd * std::abs(v.v) * v.v;
                double total_mass = v.mass + rho * v.area * v.h;

                double a = (Fw + Fb + Fd) / total_mass;
                double dhdt = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y - v.h);

                v.updateState(a, dhdt);

                currentResults.push_back(v.y);
                currentResults.push_back(v.v);
                currentResults.push_back(v.h);
            }
            results.push_back(currentResults);
        }
    }
};

#endif
