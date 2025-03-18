// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#include "ensemble.h"
#include "solver.h"
#include "interpolator.h"


// Default setup simulates a scenario of exponential tissue growth starting from a single cell.
// On a laptop with 4 cores @ 2.80GHz the default setup takes about 1 minute to run.
// There are two chemical species with concentrations c0 and c1.
// The former is produced by the starting cell, the latter at the domain boundary. 
// Both concentrations degrade linearly with degradation rates k0 and k1 respectively.
// Cells that experience a local concentration c0 below a certain threshold stop growing and differentiate. 
// The cells also slowly move towards the gradient of c1. 

int main(int argc, char* argv[]) {
    
    welcome();
    validate_parameters(); // checks that correct number of input parameters are set
    write_config(); // write cfg file to save parameters

    double Lx = 30;
    double Ly = 30;
    Domain domain(-Lx/2, -Ly/2, Lx/2, Ly/2); // initialize rectangular domain
    Ensemble ensemble("ensemble.off", domain); // init ensemble with input file

    // define reaction model R(c,k,t)
    // c: vector of concentrations, k: vector of kinetic coefficients, t: time
    Reaction linearDegradation = [](const std::vector<double>& c, const std::vector<double>& k, double t) {
        return std::vector<double> {
            - k[0] * c[0] + k[2], // linear degradation plus constant production of c0
            - k[1] * c[1]  // linear degradation of c1
        };
    }; 

    Solver solver(domain, dx, linearDegradation); // init solver
    Interpolator interpolator(ensemble, solver); // init interpolator 

    // choose method for exterior velocity interpolation (outside of polygons). only relevant if ADVECTION_DILUTION_EN is true
    interpolator.ext_interpolation_method = InterpolationMethod::ZERO;
    
    // specify custom boundary conditions for species 1. by default everything is zero-flux
    solver.boundary[1].west = {BoundaryCondition::Type::Dirichlet, 1};
    solver.boundary[1].east = {BoundaryCondition::Type::Dirichlet, 0};

    // user-defined lambda to adjust kinetic coefficients based on cell properties (e.g. cell type, position, etc). 
    ensemble.kinCoeff = [](const Polygon& self) { 
        // only starting cell (index 0) produces c0 with rate k2, make zero for other cells. leave k0/k1 as is
        return std::vector<double> {self.k[0], self.k[1], (self.polygon_index() == 0) * self.k[2]}; 
    };

    // user-defined lambdas f(c,∇c,t) for concentration effects on cell behavior
    // c: vector of concentrations, grad_c: vector of concentration gradients, t: time
    ensemble.acceleration = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        return 3e2 * grad_c[1]; // move towards gradient of c1
    };
    constexpr int T = Ns*Nf*dt; // final time
    ensemble.cellType = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (c[0] < 0.01 && t > T/2) return 1; // differentiate cell type if concentration c0 falls below threshold after T/2
        else if (self.cell_type == 1) return 1; // keep cell type once differentiated
        else return 0; 
    };
    ensemble.growthRate = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (self.cell_type == 1) return 0.0; // stop growth if cell has differentiated
        else return self.alpha; // otherwise use the default growth rate
    };
    // ensemble.maxArea =
    // ensemble.areaStiffness =
    // ensemble.lineTension = 
    // ensemble.edgeContractilityStiffness = 
    // ensemble.bendingStiffness =
    // ensemble.adhesionStiffness =
    // ensemble.dynamicFriction =
    
    // run simulation
    ensemble.output(0); // print the initial state (polygons)
    solver.output(0); // print the initial state (grid)
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(dt); // advance mechanical ensemble
            interpolator.scatter(); // interpolate data to grid
            solver.step(dt); // advance chemical solver
            interpolator.gather(); // interpolate data back to ensemble
        } 
        ensemble.output(f); // print frame
        solver.output(f); // print frame
    }
    ensemble.write_OFF("final_state.off"); // save final state, can reuse as input later
    goodbye();
}