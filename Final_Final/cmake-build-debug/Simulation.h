// Simulation.h
#ifndef SIMULATION_H
#define SIMULATION_H

#include <vector>
#include "Vessel.h"

class Simulation {
protected:
    double T, dt, g, rho, Cd; // Simulation parameters
    std::vector<std::vector<double>> results; // Store simulation results

public:
    virtual void runSimulation(std::vector<Vessel> &vessels) = 0;

    void setParameters(double totalTime, double timeStep, double gravity, double fluidDensity, double dragCoefficient) {
        T = totalTime;
        dt = timeStep;
        g = gravity;
        rho = fluidDensity;
        Cd = dragCoefficient;
    }

    const std::vector<std::vector<double>> &getResults() const {
        return results;
    }

    virtual ~Simulation() = default;
};

#endif