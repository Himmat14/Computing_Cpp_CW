// RK4Simulation.h
#ifndef RK4SIMULATION_H
#define RK4SIMULATION_H

#include "Simulation.h"

class RK4Simulation : public Simulation {
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

                double Fw1 = v.mass * g + rho * g * v.area * v.h;
                double Fb1 = -rho * g * v.area * (v.h0 + v.y - v.h);
                double Fd1 = -0.5 * rho * v.area * Cd * std::abs(v.v) * v.v;
                double total_mass1 = v.mass + rho * v.area * v.h;

                double a1 = (Fw1 + Fb1 + Fd1) / total_mass1;
                double dhdt1 = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y - v.h);

                double k1_v = a1;
                double k1_y = v.v;
                double k1_h = dhdt1;

                double Fw2 = v.mass * g + rho * g * v.area * (v.h + 0.5 * k1_h * dt);
                double Fb2 = -rho * g * v.area * (v.h0 + v.y + 0.5 * k1_y * dt - (v.h + 0.5 * k1_h * dt));
                double Fd2 = -0.5 * rho * v.area * Cd * std::abs(v.v + 0.5 * k1_v * dt) * (v.v + 0.5 * k1_v * dt);
                double total_mass2 = v.mass + rho * v.area * (v.h + 0.5 * k1_h * dt);

                double a2 = (Fw2 + Fb2 + Fd2) / total_mass2;
                double dhdt2 = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y + 0.5 * k1_y * dt - (v.h + 0.5 * k1_h * dt));

                double k2_v = a2;
                double k2_y = v.v + 0.5 * k1_v * dt;
                double k2_h = dhdt2;

                double Fw3 = v.mass * g + rho * g * v.area * (v.h + 0.5 * k2_h * dt);
                double Fb3 = -rho * g * v.area * (v.h0 + v.y + 0.5 * k2_y * dt - (v.h + 0.5 * k2_h * dt));
                double Fd3 = -0.5 * rho * v.area * Cd * std::abs(v.v + 0.5 * k2_v * dt) * (v.v + 0.5 * k2_v * dt);
                double total_mass3 = v.mass + rho * v.area * (v.h + 0.5 * k2_h * dt);

                double a3 = (Fw3 + Fb3 + Fd3) / total_mass3;
                double dhdt3 = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y + 0.5 * k2_y * dt - (v.h + 0.5 * k2_h * dt));

                double k3_v = a3;
                double k3_y = v.v + 0.5 * k2_v * dt;
                double k3_h = dhdt3;

                double Fw4 = v.mass * g + rho * g * v.area * (v.h + k3_h * dt);
                double Fb4 = -rho * g * v.area * (v.h0 + v.y + k3_y * dt - (v.h + k3_h * dt));
                double Fd4 = -0.5 * rho * v.area * Cd * std::abs(v.v + k3_v * dt) * (v.v + k3_v * dt);
                double total_mass4 = v.mass + rho * v.area * (v.h + k3_h * dt);

                double a4 = (Fw4 + Fb4 + Fd4) / total_mass4;
                double dhdt4 = (v.flow_rate_coeff / (rho * v.area)) * (v.h0 + v.y + k3_y * dt - (v.h + k3_h * dt));

                double k4_v = a4;
                double k4_y = v.v + k3_v * dt;
                double k4_h = dhdt4;

                double final_v = v.v + (dt / 6) * (k1_v + 2 * k2_v + 2 * k3_v + k4_v);
                double final_y = v.y + (dt / 6) * (k1_y + 2 * k2_y + 2 * k3_y + k4_y);
                double final_h = v.h + (dt / 6) * (k1_h + 2 * k2_h + 2 * k3_h + k4_h);

                v.v = final_v;
                v.y = final_y;
                v.h = final_h;

                currentResults.push_back(v.y);
                currentResults.push_back(v.v);
                currentResults.push_back(v.h);
            }
            results.push_back(currentResults);
        }
    }
};

#endif
