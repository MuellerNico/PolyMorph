#include <vector>
#include <functional>
#include <fstream>
#include <cmath>

// size of the FD grid. 
// x always corresponds to i and y to j.
constexpr int Nx = 3; 
constexpr int Ny = 1000;

struct Grid {
    std::vector<std::vector<double>> data;
    Grid () : data(Nx, std::vector<double>(Ny, 0.0)) {} // constructor always creates Nx x Ny grid
    double& operator()(int i, int j) { return data[i][j]; }
};

struct Reaction {
    double operator()(double u, int i, int j) {
        return 0.0; // ToDo: implement reaction function
    }
};

struct LinearDegradation : Reaction {
    double k;
    LinearDegradation(double k): k(k) {}
    double operator()(double u, int i, int j) {
        return -k * u;
    }
};

struct Solver { 
    double box_position_x;  // TODO: change this ugly datastruct. Maybe a dimensions struct
    double box_position_y;   
    Grid u; // concentration of the morphogen
    double D; // diffusion coefficient. Later also a grid datastructure
    double dx; // grid spacing
    double dt; // time step
    Reaction R; // reaction term

    // initialize the grid with a given initial condition
    Solver(const Grid u0, const double D = 1.0, const double dx = 0.1, 
            double dt = 1e-4, Reaction R = LinearDegradation(0.1)) {
        this->u = u0;
        this->D = D;
        this->dx = dx;
        this->dt = dt;
        this->R = R;
        box_position_x = - Nx/2 * dx; // midpoint at 0
        box_position_y = - Ny/2 * dx;
        enforce_bc(u);
    }

    void step() { 
        Grid unew;
        // Forward Euler with central differences ToDo: adapt for variable diffusion coefficient
        #pragma omp parallel for
        for (int i = 1; i < Nx - 1; i++) {
            for (int j = 1; j < Ny - 1; j++) {
                unew(i, j) = u(i, j) + dt * (
                    D / (dx*dx) * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1) - 4 * u(i, j))
                    + R(u(i, j), i, j)
                );
            }
        }
        enforce_bc(unew);
        u = unew; 
    }

    void enforce_bc(Grid& u) {
        // Zero-flux boundary conditions 
        for (int i = 0; i < Nx; i++) {
            u(i, 0) = 1;
            u(i, Ny - 1) = 2 * u(i, Ny - 2);
        }
        for (int j = 0; j < Ny; j++) {
            u(0, j) = 2 * u(1, j);
            u(Nx - 1, j) = 2 * u(Nx - 2, j);
        }
    }

    void output(const std::size_t frame) {    // f: frame number
        char filename [19]; 
        snprintf(filename, 19, "rd_frame%06zu.vts", frame);        
        std::ofstream file(filename);

        file << "<?xml version=\"1.0\"?>" << std::endl;
        file << "<VTKFile type=\"StructuredGrid\" version=\"0.1\">" << std::endl;
        file << "<StructuredGrid WholeExtent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Piece Extent=\"0 " << Nx-1 << " 0 " << Ny-1 << " 0 0\">" << std::endl;
        file << "<Points>" << std::endl;
        file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                double x = box_position_x + i * dx;
                double y = box_position_y + j * dx;
                file << x << " " << y << " 0" << std::endl;
            }
        }
        file << "</DataArray>" << std::endl;
        file << "</Points>" << std::endl;
        file << "<PointData Scalars=\"scalars\">" << std::endl;
        file << "<DataArray type=\"Float64\" Name=\"u\" format=\"ascii\">" << std::endl;
        for (int i = 0; i < Nx; i++) {
            for (int j = 0; j < Ny; j++) {
                file << u(i, j) << " ";
            }
            file << std::endl;
        }
        file << "</DataArray>" << std::endl;
        file << "</PointData>" << std::endl;
        file << "</Piece>" << std::endl;
        file << "</StructuredGrid>" << std::endl;
        file << "</VTKFile>" << std::endl;
        file.close();         
    }
};

// helper function for IC
/*
Grid create_gaussian() {
    Grid u;
    for(int i = 0; i < Nx; i++)
        for(int j = 0; j < Ny; j++)
            u(i, j) = std::exp(-((i - N/2)*(i - N/2) + (j - N/2)*(j - N/2)) / 100.0);
    u(N/2, N/2) = 1.0;
    return u;
}
*/