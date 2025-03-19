// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef PARAM_H
#define PARAM_H

#include <cmath>
#include <vector>

#include "geometry.h"

// -- Original PolyHoop parameters --
constexpr double h = 0.01; // [L] edge thickness
constexpr double lmin = 0.02; // [L] minimum edge length
constexpr double lmax = 0.2; // [L] maximum edge length
constexpr double C = 1; // [-] circularity index
constexpr double alpha_mu = 1; // [L^2/T] mean area growth rate
constexpr double alpha_CV = 0.3; // [-] coefficient of variation of area growth rate
constexpr double beta = 0.8; // [-] minimum area fraction for growth
constexpr double Amin = 0; // [L^2] minimum area
constexpr double Amax_mu = M_PI; // [L^2] mean maximum area
constexpr double Amax_CV = 0.1; // [-] coefficient of variation of maximum area
constexpr double gam = 1e4; // [L/T^2] line tension per vertex mass
constexpr double ka = 1e5; // [1/(L^2*T^2)] area stiffness per vertex mass
constexpr double kl = 1e4; // [L/T^2] edge contractility stiffness per vertex mass
constexpr double kb = 0; // [L^3/T^2] bending stiffness per vertex mass
constexpr double kr = 1e7; // [1/T^2] repulsion stiffness per vertex mass
constexpr double kh = 1e6; // [1/T^2] adhesion stiffness per vertex mass
constexpr double sh = 0.01; // [L] adhesion hardening zone size
constexpr double ss = 0.01; // [L] adhesion softening zone size
constexpr double mu = 0; // [-] dynamic friction coefficient
constexpr double rho = 0; // [1/L^2] fluid mass density per vertex mass
constexpr double g = 0; // [L/T^2] gravitational acceleration
constexpr double gl = 0; // [L/T^2] edge gravitational acceleration
constexpr double cv = 20; // [1/T] viscous damping rate
constexpr double cd = 0; // [-] drag coefficient
constexpr double cc = 30; // [1/T] collision damping rate
constexpr double dt = 1e-4; // [T] time step

constexpr std::size_t Nf = 100; // number of output frames
constexpr std::size_t Ns = 1000; // number of time steps between frames
constexpr std::size_t Nr = 0; // number of rigid polygons

constexpr double drmax = h + sh + ss; // maximum interaction distance

// -- New PolyMorph parameters --
constexpr std::size_t NUM_SPECIES = 2; // number of diffusable species (defines size of vectors D,c,grad_c)
constexpr std::size_t NUM_KINETIC = 3; // number of kinetic coefficients (defines size of vector k)
constexpr bool ADVECTION_DILUTION = true; // enable advection-dilution and calculate velocity field
constexpr int RNG_SEED = 0; // random number generator seed

const std::vector<double> k0 =   {1, 0.1, 0}; // background reaction coefficients
const std::vector<double> k_mu = {1, 0.1, 1}; // mean reaction coefficients
const std::vector<double> k_CV = {0.3, 0, 0}; // [-] CV of reaction coefficients
const std::vector<double> D0 =   {32, 32}; // [L^2/T] background diffusivity (recommended to not be zero)
const std::vector<double> D_mu = {32, 32}; // [L^2/T] mean diffusivities
const std::vector<double> D_CV = {0.3, 0.3}; // [-] CV of diffusivities
const std::vector<double> anisotropy = {1, 1}; // [-] diffusion anisotropy (default 1 = isotropic)

Domain domain = {-15, -15, 15, 15}; // rectangular spatial domain (x0, y0, x1, y1)
constexpr double dx = 0.3; // [L] grid spacing for solver (same in x and y directions)
constexpr double dist_cutoff = 2; // [-] used to cut off lognormal distributions at mu*cutoff (only used for diffusivities currently)
constexpr double kd = kr; // [1/T^2] domain boundary stiffness (too high can cause instabilities when tissue fills out entire domain)
const double IDW_rcut = 4 * std::sqrt(Amax_mu); // [L] cutoff radius for external IDW interpolation of velocity field

std::string input_file = "ensemble.off"; // default input ensemble file
std::string output_folder = "output"; // default output folder

// choose what to write to VTK files for visualization
namespace Output {
    constexpr bool c = true; // concentration
    constexpr bool grad_c = true; // concentration gradient
    constexpr bool D = true; // diffusion coefficient
    constexpr bool k = true; // kinetic coefficients
    constexpr bool parent_idx = true; // polygon (cell) index
    constexpr bool cell_type = true; // cell type
    constexpr bool v = true; // velocity field
};

#endif // PARAM_H
