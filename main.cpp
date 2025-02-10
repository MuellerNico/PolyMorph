#include "ensembleController.h"

int main(int argc, char* argv[]) {
    welcome();
    validate_parameters(); // checks that correct number of kinetic parameters are set
    write_config(); // write cfg file to save parameters
    
    double L = 30;
    Domain domain(-L/2, -L/2, L/2, L/2); // initialize rectangular domain
    Ensemble ensemble("ensemble.off", domain); // init ensemble with input file
    Reaction reaction = Reactions::linearDegradation; // define reaction model
    Solver solver(domain, dx, reaction); // init solver
    Interpolator interpolator(ensemble, solver); // init interpolator
    interpolator.ext_interpolation_method = InterpolationMethod::ZERO; // interpolation for exterior nodes (interior is always IDW)
    
    // set boundary conditions (default: zero-flux)
    //solver.boundary[1].west = {BoundaryCondition::Type::Dirichlet, 0.01}; 
    // define production lambda (default no cell-based production)
    ensemble.is_producing = [](const Polygon& p) { 
       //return std::vector<bool> {p.global_index() == 0, false, false};
       return std::vector<bool>(NUM_SPECIES, p.global_index() == 0); // starting cell (index 0) produces
    }; 
    // define lambdas for concentration effects on cell behavior
    ensemble.cellTypeEffect = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (c[0] < 0.005) return 1; // differentiate cell type if concentration falls below threshold
        //else if (c[0] + c[1] + c[2] > 0.6) return 2;
        else return 0; 
    };
    ensemble.growthRateEffect = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        if (c[0] < 0.005) return 0.0; // stop growth if concentration falls below threshold
        else return self.alpha; 
    };
    ensemble.accelerationEffect = [](const Polygon& self, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { 
        return Point(0, 0); // no acceleration effect
    };
    
    ensemble.output(0); // print the initial state
    solver.output(0); // print the initial state
    for (std::size_t f = 1; f <= Nf; ++f) {
        for (std::size_t s = 0; s < Ns; ++s) {
            ensemble.step(); // advance mechanical ensemble
            interpolator.scatter(); // interpolate data to grid
            solver.step(dt); // advance chemical solver
            interpolator.gather(); // interpolate data back to ensemble
            //domain.step(dt); // only needed if domain is growing/shrinking
        } 
        ensemble.output(f); // print frame
        solver.output(f); // print frame
    }
    //ensemble.write_OFF("final_state.off"); // save final state, can reuse as input later
    goodbye();
}