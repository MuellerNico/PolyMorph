#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <functional>
#include <fstream>
#include <iostream>

#include "grid.h"
#include "geometry.h"
#include "domain.h"

using Reaction = std::function<std::vector<double>(const std::vector<double>& c, const std::vector<double>& k, double t)>;

struct BoundaryCondition {
    enum class Type {
        Dirichlet,
        Neumann
    };

    Type type; // Dirichlet or Neumann
    double value; // value at boundary (Dirichlet case) or derivative at boundary (Neumann case)
};

struct Boundary { // boundary conditions for one species
    BoundaryCondition north, south, east, west;
    
    static Boundary zeroFlux() { // default boundary condition
        return {
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0},
            {BoundaryCondition::Type::Neumann, 0}
        };
    }

    static std::vector<Boundary> zeroFlux(size_t n) { // default boundary conditions for n species
        return std::vector<Boundary>(n, zeroFlux());
    }
};

struct Solver { 
    Domain& domain; // rectangular simulation domain
    int Nx, Ny; // number of grid points
    double dx; // grid spacing
    double t; // time
    std::vector<Boundary> boundary; // boundary conditions (one Boundary per species)
    Reaction R; // reaction term
    Grid<int> parent_idx; // polygon idx
    Grid<std::vector<double>> c; // concentrations
    Grid<std::vector<double>> cnew; // temporary grid for updating concentrations
    Grid<std::vector<double>> D; // diffusion coefficients
    Grid<std::vector<double>> k; // kinetic coefficients
    Grid<Point> velocity; // velocity field
    Grid<std::vector<Point>> grad_c; // concentration gradient

    Solver(Domain& domain, const double dx, Reaction R) : t(0), domain(domain), R(R), dx(dx) {
        this->boundary = Boundary::zeroFlux(NUM_SPECIES);
        this->Nx = domain.width() / dx + 1;
        this->Ny = domain.height() / dx + 1;
        this->c = Grid<std::vector<double>>(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));
        this->cnew = Grid<std::vector<double>>(Nx, Ny, std::vector<double>(NUM_SPECIES, 0.0));
        this->grad_c = Grid<std::vector<Point>>(Nx, Ny, std::vector<Point>(NUM_SPECIES, Point(0, 0)));
        this->R = [](const std::vector<double>& c, const std::vector<double>& k, double t) { return std::vector<double>(NUM_SPECIES, 0.0); };
        
        std::cout << "solver dimensions: Nx=" << Nx << " Ny=" << Ny << " dx=" << dx << std::endl;
        std::cout << "domain: [" << domain.x0 << ", " << domain.x1 << "] x [" << domain.y0 << ", " << domain.y1 << "]" << std::endl;

        // initialize with background values
        parent_idx = Grid<int>(Nx, Ny, -2);
        D = Grid<std::vector<double>>(Nx, Ny, D0);
        k = Grid<std::vector<double>>(Nx, Ny, k0);
        
        if (ADVECTION_DILUTION_EN) {
            velocity = Grid<Point>(Nx, Ny, Point(0, 0));
        }
    }

    /// @param dt time step: allows for independent time stepping between ensemble and solver
    void step(double dt) {
        // precompute
        double two_dx = 2 * dx; 
        double dx2 = dx * dx; 
        // Forward Euler with central differences
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {   
                const std::vector<double> reaction = R(c(i, j), k(i, j), t);
                for (int sp = 0; sp < NUM_SPECIES; sp++) {
                    // dirichlet boundary conditions
                    if      (i == 0     && boundary[sp].west.type  == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].west.value; continue; } 
                    else if (i == Nx-1  && boundary[sp].east.type  == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].east.value; continue; } 
                    else if (j == 0     && boundary[sp].south.type == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].south.value; continue; } 
                    else if (j == Ny-1  && boundary[sp].north.type == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].north.value; continue; } 
                    else {
                        // account for neumann BDC
                        const double n = (j == Ny-1) ? c(i, j-1)[sp] + two_dx*boundary[sp].north.value : c(i, j+1)[sp]; 
                        const double s = (j == 0)    ? c(i, j+1)[sp] - two_dx*boundary[sp].south.value : c(i, j-1)[sp];
                        const double e = (i == Nx-1) ? c(i-1, j)[sp] + two_dx*boundary[sp].east.value  : c(i+1, j)[sp];
                        const double w = (i == 0)    ? c(i+1, j)[sp] - two_dx*boundary[sp].west.value  : c(i-1, j)[sp];
                        // calculate diffusion term
                        const double diffusion = D(i, j)[sp] / dx2 * (e + w + anisotropy[sp] * (n + s) - 2 * (1 + anisotropy[sp]) * c(i, j)[sp]);
                        // update grid point
                        cnew(i, j)[sp] = c(i, j)[sp] + dt * (diffusion + reaction[sp]); 
                        grad_c(i, j)[sp] = {(e - w) / two_dx, (n - s) / two_dx};
                        
                        if(ADVECTION_DILUTION_EN) {
                            const double advection = velocity(i, j) * grad_c(i, j)[sp]; // dot product
                            const double dvdx = (j == Ny-1 || j == 0) ? 0 : (velocity(i, j+1).y - velocity(i, j-1).y) / two_dx;
                            const double dvdy = (i == Nx-1 || i == 0) ? 0 : (velocity(i+1, j).x - velocity(i-1, j).x) / two_dx;
                            const double dilution = c(i, j)[sp] * (dvdx + dvdy);
                            // update grid point
                            cnew(i, j)[sp] -= dt * (advection + dilution);
                        }                    
                    }
                }
            }
        }
        // update state
        c.parallel_copy_from(cnew); // parallelized assignment operator
        t += dt; // advance the time
    }

    void output(const std::size_t frame) {
        char filename [19]; 
        snprintf(filename, 19, "rd_frame%06zu.vts", frame);        
        std::ofstream file(filename);

        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
        file << "<StructuredGrid WholeExtent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Piece Extent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        // define points
        file << "<Points>" << std::endl;
        file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = domain.x0 + i * dx;
                double y = domain.y0 + j * dx;
                file << x << " " << y << " 0" << std::endl;
            }
        }
        file << "</DataArray>" << std::endl;
        file << "</Points>" << std::endl;
        file << "<PointData Scalars=\"scalars\">" << std::endl; // start point data
        
        if (Output::c) file << c.to_vtk("c");
        if (Output::grad_c) file << grad_c.to_vtk("gradient");
        if (Output::parent_idx) file << parent_idx.to_vtk("parent_idx");
        if (Output::D) file << D.to_vtk("D");
        if (Output::k) file << k.to_vtk("k");
        if (Output::velocity && ADVECTION_DILUTION_EN) file << velocity.to_vtk("velocity");
        
        file << "</PointData>" << std::endl;    // end of point data
        file << "</Piece>" << std::endl;
        file << "</StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();          
    }
};

#endif