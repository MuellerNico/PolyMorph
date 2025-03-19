// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#include "ensemble.h"
#include "solver.h"
#include "interpolator.h"

// The default setup simulates a scenario of exponential tissue growth starting from a single cell.
// On a laptop with 4 cores at 2.80GHz, the default setup takes about 1 minute to run.
// There are two chemical species with concentrations c1 and c2.
// The concentrations degrade linearly with degradation rates k1 and k2 respectively.
// The first is produced with rate k3 by the starting cell, the second at the left domain boundary.
// Cells that experience a local concentration c1 below a certain threshold after some time differentiate and stop growing.
// The cells migrate up the gradient of c2.

int main(int argc, char* argv[]) {
    // get input file from optional command line argument
    if (argc > 1) {
        input_file = argv[1];
    }
    
    // get output folder from optional command line argument
    if (argc > 2) {
        output_folder = argv[2];
    }
    
    welcome();
    check_parameters(); // check that correct number of input parameters are set
    
    std::system((std::string("mkdir -p ") + output_folder).c_str()); // create the output folder if it doesn't exist yet
    
    write_config(); // write cfg file to save parameters for reproducibility
    
    Ensemble ensemble(input_file.c_str()); // initialize ensemble

    // define reaction model R(p,c,k,t)
    // p: polygon index, c: vector of concentrations, k: vector of kinetic coefficients, t: time
    Reaction R = [&ensemble](std::size_t p, const std::vector<double>& c, const std::vector<double>& k, double t) {
        return std::vector<double> {
            - k[0] * c[0] + (ensemble.polygons[p].cell_type == 0) * k[2], // linear degradation plus constant production of c1 in the starting cell
            - k[1] * c[1]  // linear degradation of c2
        };
    };

    Solver solver(R); // initialize solver
    Interpolator interpolator(ensemble, solver); // initialize interpolator
    
    // if enabling advection-dilution, choose method for exterior velocity interpolation
    interpolator.exterior_method = InterpolationMethod::ZERO; // one of IDW, BILINEAR, ZERO

    // specify custom boundary conditions for species 2, everything else is zero-flux
    solver.boundary[1].west = {BoundaryCondition::Type::Dirichlet, 1};

    // user-defined lambdas f(c,∇c,t) for concentration effects on cell behavior
    // p: polygon, c: vector of concentrations, grad_c: vector of concentration gradients, t: time
    ensemble.acceleration = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) {
        return 5e2 * grad_c[1]; // chemotaxis upward the gradient of c2
    };
    constexpr int T = Ns*Nf*dt; // final time
    ensemble.cellType = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) {
        if (p.index() == 0) return 0; // first cell remains special
        if (c[0] < 0.01 && t > T/2) return 2; // differentiate if concentration c1 falls below a threshold after T/2
        else if (p.cell_type == 2) return 2; // once differentiated, keep cell type
        else return 1; // growing
    };
    ensemble.growthRate = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) {
        if (p.cell_type == 2) return 0.; // differentiated cells are quiescent
        else return p.alpha; // all others use the default growth rate
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
        ensemble.output(f); // print frame (polygons)
        solver.output(f); // print frame (grid)
    }
    ensemble.write_OFF("final_state.off"); // save final state, can be reused as input later
    goodbye();
}
