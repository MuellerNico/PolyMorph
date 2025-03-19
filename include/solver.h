// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <functional>
#include <fstream>
#include <iostream>

#include "grid.h"
#include "geometry.h"
#include "param.h"

using Reaction = std::function<std::vector<double>(std::size_t p, const std::vector<double>& c, const std::vector<double>& k, double t)>;

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

    static std::vector<Boundary> zeroFlux(std::size_t n) { // default boundary conditions for n species
        return std::vector<Boundary>(n, zeroFlux());
    }
};

struct Solver {
    int Nx, Ny; // number of grid points
    double t; // time
    std::vector<Boundary> boundary; // boundary conditions (one Boundary object per species)
    Reaction R; // reaction term
    Grid<int> parent_idx; // polygon index
    Grid<std::vector<double>> c; // concentrations
    Grid<std::vector<double>> cnew; // temporary grid for updating concentrations
    Grid<std::vector<double>> D; // diffusion coefficients
    Grid<std::vector<double>> k; // kinetic coefficients
    Grid<Point> v; // velocity field
    Grid<std::vector<Point>> grad_c; // concentration gradient

    Solver(Reaction R) : t(0), R(R) {
        boundary = Boundary::zeroFlux(NUM_SPECIES);
        Nx = (domain.x1 - domain.x0) / dx + 1;
        Ny = (domain.y1 - domain.y0) / dx + 1;
        c = Grid<std::vector<double>>(Nx, Ny, std::vector<double>(NUM_SPECIES, 0));
        cnew = Grid<std::vector<double>>(Nx, Ny, std::vector<double>(NUM_SPECIES, 0));
        grad_c = Grid<std::vector<Point>>(Nx, Ny, std::vector<Point>(NUM_SPECIES, {0, 0}));

        // initialize with background values
        parent_idx = Grid<int>(Nx, Ny, -2);
        D = Grid<std::vector<double>>(Nx, Ny, D0);
        k = Grid<std::vector<double>>(Nx, Ny, k0);

        if (ADVECTION_DILUTION) {
            v = Grid<Point>(Nx, Ny, {0, 0}); // init to zero (also serves as boundary condition)
        }
    }

    void step(double dt) {
        static constexpr double two_dx = 2 * dx;
        // Forward Euler with central differences
        #pragma omp parallel for collapse(2)
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                const std::vector<double> reaction = R(parent_idx(i, j), c(i, j), k(i, j), t);
                for (std::size_t sp = 0; sp < NUM_SPECIES; sp++) {
                    // Dirichlet boundary conditions
                    if      (i == 0     && boundary[sp].west.type  == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].west.value; continue; }
                    else if (i == Nx-1  && boundary[sp].east.type  == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].east.value; continue; }
                    else if (j == 0     && boundary[sp].south.type == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].south.value; continue; }
                    else if (j == Ny-1  && boundary[sp].north.type == BoundaryCondition::Type::Dirichlet) { cnew(i, j)[sp] = boundary[sp].north.value; continue; }
                    else {
                        // account for Neumann boundary conditions
                        const double n = (j == Ny-1) ? c(i, j-1)[sp] + two_dx * boundary[sp].north.value : c(i, j+1)[sp];
                        const double s = (j == 0)    ? c(i, j+1)[sp] - two_dx * boundary[sp].south.value : c(i, j-1)[sp];
                        const double e = (i == Nx-1) ? c(i-1, j)[sp] + two_dx * boundary[sp].east.value  : c(i+1, j)[sp];
                        const double w = (i == 0)    ? c(i+1, j)[sp] - two_dx * boundary[sp].west.value  : c(i-1, j)[sp];
                        // calculate diffusion term
                        const double diffusion = D(i, j)[sp] / (dx * dx) * (e + w + anisotropy[sp] * (n + s) - 2 * (1 + anisotropy[sp]) * c(i, j)[sp]);
                        // update grid point
                        cnew(i, j)[sp] = c(i, j)[sp] + dt * (diffusion + reaction[sp]);
                        grad_c(i, j)[sp] = {(e - w) / two_dx, (n - s) / two_dx};

                        if (ADVECTION_DILUTION) {
                            const double advection = v(i, j) * grad_c(i, j)[sp]; // dot product
                            const double dvdx = (j == Ny-1 || j == 0) ? 0 : (v(i, j+1).y - v(i, j-1).y) / two_dx;
                            const double dvdy = (i == Nx-1 || i == 0) ? 0 : (v(i+1, j).x - v(i-1, j).x) / two_dx;
                            const double dilution = c(i, j)[sp] * (dvdx + dvdy);
                            cnew(i, j)[sp] -= dt * (advection + dilution); // update grid point
                        }
                    }
                }
            }
        }
        c.swap(cnew); // update the state
        t += dt; // advance the time
    }

    void output(const std::size_t frame) {
        char filename [21];
        snprintf(filename, 21, "grid_frame%06zu.vts", frame);
        std::ofstream file(output_folder + "/" + filename);
        if (!file) {
            std::cerr << "Error: Could not open file " << output_folder + "/" << filename << std::endl;
            return;
        }
        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
        file << "    <StructuredGrid WholeExtent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "        <Piece Extent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "            <Points>" << std::endl;
        file << "                <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                const double x = domain.x0 + i * dx;
                const double y = domain.y0 + j * dx;
                file << x << " " << y << " 0" << std::endl;
            }
        }
        file << "                </DataArray>" << std::endl;
        file << "            </Points>" << std::endl;
        file << "            <PointData>" << std::endl;
        if (Output::c) file << c.to_vtk("c");
        if (Output::grad_c) file << grad_c.to_vtk("grad_c");
        if (Output::parent_idx) file << parent_idx.to_vtk("parent_idx");
        if (Output::D) file << D.to_vtk("D");
        if (Output::k) file << k.to_vtk("k");
        if (Output::v && ADVECTION_DILUTION) file << v.to_vtk("v");
        file << "            </PointData>" << std::endl;
        file << "        </Piece>" << std::endl;
        file << "    </StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();
    }
};

#endif
