// Cpp CW
//Himmat Kaul 02376386

//include packages and libraries
#include <vector>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// define vessel class
class Vessel {
public:
    //define vessel parameters
    double mass;
    double area;
    double flow_rate_coeff;
    double y0;
    double v0;
    double h0;
    //initialise dynamic response properties
    double y;
    double v;
    double h;

    //define time step
    double delta_time;

    // Define function to update the values
    void new_vals(double accel, double WaterLevel_rateChange) {
    v += accel * delta_time;
    y += v * delta_time;
    h += WaterLevel_rateChange * delta_time;
    }

    // default constructor for vesssel classs
    Vessel(double m, double a, double k, double y_init, double v_init, double h_init)
            : mass(m), area(a), flow_rate_coeff(k), y0(y_init), v0(v_init), h0(h_init), y(y_init), v(v_init), h(h_init), delta_time(0) {}

};

//define sim parent class
class sim {
public:

    // starts the simulation takes in the vessel class
    virtual void start_sim(std::vector<class Vessel> &vessels) = 0;

    // Initialises sim global parameters
    void Initialise_parameters(double time, double Delta_Time, double g, double Density, double drag_coeff);

    // gets output from sim
    std::vector<std::vector<double>> getResults();

    //defines some protected atributes ie the global params
protected:
    double time_period;
    double delta_time;
    double g_const;
    double density;
    double Cd;

    // Store simulation results
    std::vector<std::vector<double>> results;
};

//Function to define parameters for sim
void sim::Initialise_parameters(double time, double Delta_Time, double g, double Density, double drag_coeff) {
    time_period = time;
    delta_time = Delta_Time;
    g_const = g;
    density = Density;
    Cd = drag_coeff;
}

//function defenition outside of the class for tidyness
std::vector<std::vector<double>> sim::getResults() {
    return results;
}

// Derived Class for Explicit Euler_sim sim__ child class of sim class
class Euler_sim : public sim {
//has acces to public atributes from sim

//defines protected atributes
protected:

    double Euler_intermediate(double vmass, double varea, double vh, double vh0, double vy, double vv) {

        //calucualtes forces from handout
        double Fwn = (vmass + density * varea * vh ) * g_const;
        double Fbn = -1 * (density) * (g_const) * ( varea * vh0 + varea * vy - varea * vh );
        double Fdn = varea * (-0.5) * density * Cd * std::abs(vv) * vv;

        // uses Newtons 2nd law to find the acceleration
        double a = (Fwn + Fbn + Fdn) / (vmass + density * varea * vh);
        //returns a
        return a;
    }

public:
    void start_sim(std::vector<Vessel> &vessels);
};

// defines start sim function for euler outside of class
void Euler_sim::start_sim(std::vector<Vessel> &vessels) {

    // loops for all time steps
    for (int i = 0; i <= (vessels.size() - 1); ++i) {
        vessels[i].delta_time = delta_time;
    }

    // predefines n for itteration
    double n = time_period / delta_time;

    // loops for all time steps
    for (int step = 0; step < (n + 1); ++step) {
        //defines the curreent results vector
        std::vector<double> currentResults;

        // calls for the vessel parameters one by one by passby referencing
        for (auto &v: vessels) {


            //finds acceleration using intermediate function
            double a = Euler_intermediate(v.mass, v.area, v.h, v.h0, v.y, v.v);

            // finds water level rate change through formula from handout
            double dhdt = (v.flow_rate_coeff / (density * v.area)) * (v.h0 + v.y - v.h);

            //updates the current values for the sim
            v.new_vals(a, dhdt);

            //adds the results to current results which is later stored in output file
            currentResults.push_back(v.y);
            currentResults.push_back(v.v);
            currentResults.push_back(v.h);
        }
        // adds the current results to output results vector
        results.push_back(currentResults);
    }
}

// Derived Class for rk4
class RK4sim : public sim {

    // defines public attributes
public:
    //defines public function to calculate forces
    std::vector<double> RK4_calculate_forces(double vmass, double varea, double vh, double vh0, double vy, double vv, double delta_time, double kn_v, double kn_y, double kn_h) {

        //calcualtes forves using the fomunals from handout
        double Fwn = vmass * g_const + density * g_const * varea * (vh + 0.5 * kn_h * delta_time);
        double Fbn = -density * g_const * varea * (vh0 + vy + 0.5 * kn_y * delta_time - (vh + 0.5 * kn_h * delta_time));
        double Fdn = -0.5 * density * varea * Cd * std::abs(vv + 0.5 * kn_v * delta_time) * (vv + 0.5 * kn_v * delta_time);
        double total_mass_n = vmass + density * varea * (vh + 0.5 * kn_h * delta_time);

        // returns forces as a vector
        return {Fwn, Fbn, Fdn, total_mass_n};
    }

    // defines a weigthed sum function to add up RK4 weights
    std::vector<double> Weighted_sum(double vv, double vy, double vh, std::vector<double> k1, std::vector<double> k2, std::vector<double> k3, std::vector<double> k4) {

        //calculates final avearges using intermediate weights and previous values
        double final_v = vv + (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]) * (delta_time / 6);
        double final_y = vy + (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) * (delta_time / 6);
        double final_h = vh + (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) * (delta_time / 6);

        // returns the valls for v, y , h after time step
        return {final_v,final_y,final_h};
    }

    //defines sim start function
    void start_sim(std::vector<Vessel> &vessels);
};


// defines start sim as atribute of class RK4sim
void RK4sim::start_sim(std::vector<Vessel> &vessels) {

    // loops for all vessels in the input file
    for (auto &vessel: vessels) {
        vessel.delta_time = delta_time;
    }

    // initialises the number of steps
    double steps = (time_period / delta_time);

    //loops for all steps
    for (int step = 0; step < (steps + 1); ++step) {
        // defines current results
        std::vector<double> currentResults;

        // lopps for all vessles
        for (auto & v : vessels) {

            // calls calculate forces to get forces for each time step
            //calculates forces for 1st weight
            std::vector<double> inter =  RK4_calculate_forces(v.mass, v.area, v.h, v.h0, v.y, v.v, delta_time, 0, 0, 0);

            //calcluates acceleration and the rate change in water level
            double a1 = (inter[0] + inter[1] + inter[2]) / inter[3];
            double dhdt1 = (v.flow_rate_coeff / (density * v.area)) * (v.h0 + v.y - v.h);

            //calculates weights from first pass
            double k1_v = a1;
            double k1_y = v.v;
            double k1_h = dhdt1;
            //defines weights vector
            std::vector<double> k1 = {k1_v, k1_y, k1_h};

            //2nd step
            //calculates forces for 2md weight
            inter =  RK4_calculate_forces(v.mass, v.area, v.h, v.h0, v.y, v.v, delta_time, k1_v, k1_y, k1_h);

            double a2 = (inter[0] + inter[1] + inter[2]) / inter[3];

            double dhdt2 = (v.flow_rate_coeff / (density * v.area)) * (v.h0 + v.y + 0.5 * k1_y * delta_time - (v.h + 0.5 * k1_h * delta_time));

            //calculates the weights
            double k2_v = a2;
            double k2_y = v.v + 0.5 * k1_v * delta_time;
            double k2_h = dhdt2;
            //defines weights vector
            std::vector<double> k2 = {k2_v, k2_y, k2_h};

            //3rd step
            //calculates forces for 3rd weight
            inter =  RK4_calculate_forces(v.mass, v.area, v.h, v.h0, v.y, v.v, delta_time, k2_v, k2_y, k2_h);

            double a3 = (inter[0] + inter[1] + inter[2]) / inter[3];
            double dhdt3 = (v.flow_rate_coeff / (density * v.area)) * (v.h0 + v.y + 0.5 * k2_y * delta_time - (v.h + 0.5 * k2_h * delta_time));

            //calculates the weights
            double k3_v = a3;
            double k3_y = v.v + 0.5 * k2_v * delta_time;
            double k3_h = dhdt3;
            //defines weights vector
            std::vector<double> k3 = {k3_v, k3_y, k3_h};

            //4th step
            //calculates forces for 4th weight
            inter =  RK4_calculate_forces(v.mass, v.area, v.h, v.h0, v.y, v.v, delta_time, k3_v, k3_y, k3_h);

            double a4 = (inter[0] + inter[1] + inter[2]) / inter[3];
            double dhdt4 = (v.flow_rate_coeff / (density * v.area)) * (v.h0 + v.y + k3_y * delta_time - (v.h + k3_h * delta_time));

            //calculates the weights
            double k4_v = a4;
            double k4_y = v.v + k3_v * delta_time;
            double k4_h = dhdt4;
            //defines weights vector
            std::vector<double> k4 = {k4_v, k4_y, k4_h};

            //Weighted sum step
            std::vector<double> final = Weighted_sum(v.v,v.y,v.h,k1,k2,k3,k4);

            // re-initialising final values to start values
            v.v = final[0];
            v.y = final[1];
            v.h = final[2];

            // adds results  to current results
            currentResults.push_back(v.y);
            currentResults.push_back(v.v);
            currentResults.push_back(v.h);
        }
        // adds the final time step vals to the current results
        results.push_back(currentResults);
    }
}


// Read the txt file and cdistributes params ect..

void readParameters(const std::string &filename, int &scheme, double &time_period, double &delta_time, double &g_const, double &density, double &Cd, std::vector<Vessel> &vessels) {

    //using input finle stream
    std::ifstream infile(filename);

    // validates there is a file present that can be read
    if (!infile) {
        std::cerr << "File cannot be opend!" << filename << std::endl;
        exit(1);
    }

    // if file read succesfully
    else {
        // defines a line as a string
        std::string line;
        //defines bool variable checker for reading global params
        bool readGlobalParams = false;

        // loops while there is a line to read
        while (std::getline(infile, line)) {

            // if line has hash for first charactre skips line as per handout
            if (line[0] == '#')
                continue;


            std::istringstream iss(line);


            //reads global params
            if (!readGlobalParams) {
                //checks if iss foramt is corect to read params
                if (!(iss >> scheme >> time_period >> delta_time >> g_const >> density >> Cd)) {
                    std::cerr << "Invalid Format for global parameters." << std::endl;
                    exit(1);
                }
                // sets true if read succesfully
                readGlobalParams = true;

                // reads the vessel params
            } else {
                // defines vessel params
                double m, a, k, y_init, v_init, h_init;

                //checks line format
                if (!(iss >> m >> a >> k >> y_init >> v_init >> h_init)) {
                    std::cerr << "Invalid Format for vessel parameters." << std::endl;
                    exit(1);
                }
                // adds the vessel paranms to thhe vessel vector
                vessels.push_back(Vessel(m, a, k, y_init, v_init, h_init));
            }
        }
        // checks is there is global params
        if (!readGlobalParams) {
            std::cerr << "Global parameters were not found." << std::endl;
            exit(1);
        }

        //checks for vessels
        if (vessels.empty()) {
            std::cerr << "Vessel data not found." << std::endl;
            exit(1);
        }
    }
}

// functioon to create output txt file
void generate_output(const std::string &filename, const std::vector<std::vector<double>> &results, double delta_time) {
    // creates output file
    std::ofstream outfile(filename);

    if (!outfile) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        exit(1);
    }

    // adds the values to text file
    for (int i = 0; i < results.size(); ++i) {
        outfile << std::fixed;
        outfile << i * delta_time;
        for (auto &val : results[i]) {
            outfile << " " << val;
        }
        //creates new line for next time step
        outfile << "\n";
    }
}


int main() {
    int scheme = -1;
    // initialises an errorsome scheme rewrite
    //defines the global params
    double time_period, delta_time, g_const, density, Cd;
    //defines vessel vector
    std::vector<Vessel> vessels;



    // ccaaalls the reread text file
    readParameters("parameters.txt", scheme, time_period, delta_time, g_const,
                   density, Cd, vessels);

    // creates simulation objects
    Euler_sim eulerS;
    RK4sim rk4S;

    //runns for FE
    if (scheme == 0) {
        std::cout<< "Scheme = 0, \nSimulating Explicit Euler ODE solver..."<<std::endl;
        std::cout << "Global Parameters: \nT = "<<time_period<<", dt = "<<delta_time<<", g = "<<g_const<<", rho = "<<density<<", Cd = "<<Cd<<std::endl;
        eulerS.Initialise_parameters(time_period, delta_time, g_const, density, Cd);
        eulerS.start_sim(vessels);

        //calls generate output file
        generate_output("output.txt", eulerS.getResults(), delta_time);
    }

    // runns for RK4
    else if (scheme == 1) {
        std::cout<< "Scheme = 1, Simulating Runge-Kutta 4th Order ODE solver"<<std::endl;
        std::cout << "Global Parameters: \nT = "<<time_period<<", dt = "<<delta_time<<", g = "<<g_const<<", rho = "<<density<<", Cd = "<<Cd<<std::endl;
        rk4S.Initialise_parameters(time_period, delta_time, g_const, density, Cd);
        rk4S.start_sim(vessels);

        //calls generate output file
        generate_output("output.txt", rk4S.getResults(), delta_time);
    }

    // iff incorect or no scheme read with output messagee
    else {
        std::cerr << "No Scheme Selected: \n[0] Explicit Euler \n[1] Runge-Kutta 4th order" << std::endl;
    }\

    return 0;
}