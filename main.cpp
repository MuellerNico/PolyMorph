#include "ensemble.h"
#include "solver.h"
#include "interpolator.h"

/*
* Default setup simulates a scenario of exponential tissue growth starting from a single cell.
* There are two chemical species with concentrations c0 and c1.
* The former is produced by the starting cell, the latter at the boundary. 
* Both concentrations degrade linearly with degradation rates k0 = k1 = 1.
* Cells that experience a local concentration c0 below a certain threshold stop growing and differentiate. 
* The cells also slowly move towards the gradient of c1. 
*/
int main(int argc, char* argv[]) {
    welcome();
    validate_parameters(); // checks that correct number of input parameters are set
    write_config(); // write cfg file to save parameters
    
    double Lx = 30;
    double Ly = 30;
    Domain domain(-Lx/2, -Ly/2, Lx/2, Ly/2); // initialize rectangular domain
    Ensemble ensemble("ensemble.off", domain); // init ensemble with input file

    // define reaction model R(c,k,t)
    Reaction linearDegradation = [](const std::vector<double>& c, const std::vector<double>& k, double t) {
        return std::vector<double> {
            - k[0] * c[0] + k[2], // linear degradation plus constant production of c0
            - k[1] * c[1]  // linear degradation of c1
        };
    }; 

    Solver solver(domain, dx, linearDegradation); // init solver
    Interpolator interpolator(ensemble, solver); // init interpolator
    
    // specify boundary conditions for species 1. default is zero-flux.
    solver.boundary[1].west = {BoundaryCondition::Type::Dirichlet, 1};
    solver.boundary[1].east = {BoundaryCondition::Type::Dirichlet, 0};
    
    ensemble.polygons[0].k[2] = 1; // set constant production rate for starting cell only

    // user defined lambdas f(c,∇c,t) for concentration effects on cell behavior
    ensemble.accelerationEffect = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        return 3e2 * grad_c[1]; // move towards gradient of c1
    };
    int T = Ns*Nf*dt; // final time
    ensemble.cellTypeEffect = [T](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (c[0] < 0.005 && t > T/2) return 1; // differentiate cell type if concentration falls below threshold after T/2
        else if (self.cell_type == 1) return 1; // keep cell type once differentiated
        else return 0; 
    };
    ensemble.growthRateEffect = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (self.cell_type == 1) return 0.0; // stop growth if cell has differentiated
        else return self.alpha; 
    };
    
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); // advance mechanical ensemble
            interpolator.scatter(); // interpolate data to grid
            solver.step(dt); // advance chemical solver
            interpolator.gather(); // interpolate data back to ensemble
        } 
        ensemble.output(f); // print frame
        solver.output(f); // print frame
    }
    //ensemble.write_OFF("final_state.off"); // save final state, can reuse as input later
    goodbye();
}