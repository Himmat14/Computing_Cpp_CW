// Vessel.h
#ifndef VESSEL_H
#define VESSEL_H

class Vessel {
public:
    double mass, area, flow_rate_coeff, y0, v0, h0; // Initial properties
    double y, v, h; // Current state (position, velocity, water level)
    double dt; // Time step for updates

    Vessel(double m, double a, double k, double y_init, double v_init, double h_init)
            : mass(m), area(a), flow_rate_coeff(k), y0(y_init), v0(v_init), h0(h_init), y(y_init), v(v_init), h(h_init), dt(0) {}

    void updateState(double acceleration, double rateOfChangeWaterLevel) {
        v += acceleration * dt;
        y += v * dt;
        h += rateOfChangeWaterLevel * dt;
    }
};

#endif
