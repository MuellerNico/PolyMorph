// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef UTILS_H
#define UTILS_H

#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>
#include <random>
#include <omp.h>
#include <ctime>
#include <iomanip>

#include "param.h"

std::string getCurrentDateTime() {
    auto time = std::time(nullptr);
    auto loc_time = *std::localtime(&time);
    std::ostringstream oss;
    oss << std::put_time(&loc_time, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

void check_parameters() {
    assert(D_mu.size() == NUM_SPECIES && D_CV.size() == NUM_SPECIES);
    assert(D0.size() == NUM_SPECIES);
    assert(k_mu.size() == NUM_KINETIC && k_CV.size() == NUM_KINETIC);
    assert(k0.size() == NUM_KINETIC);
    assert(anisotropy.size() == NUM_SPECIES);
}

template <typename T>
std::string to_string(const std::vector<T>& v) {
    std::string s = "[";
    for (const T& x : v) {
        s += std::to_string(x) + ", ";
    }
    s.pop_back();
    s.pop_back();
    s += "]";
    return s;
}

// create a vector of lognormal distributions from vectors of means and CVs
std::vector<std::lognormal_distribution<>> create_lognormal(const std::vector<double>& mu, const std::vector<double>& CV) {
    assert(mu.size() == CV.size());
    std::vector<std::lognormal_distribution<>> dists;
    for (unsigned i = 0; i < mu.size(); i++) {
        const double sigma = std::sqrt(std::log(1 + CV[i] * CV[i]));  // Standard deviation of the underlying normal
        const double mu_prime = std::log(mu[i]) - 0.5 * sigma * sigma;  // Adjust mean of the underlying normal
        dists.push_back(std::lognormal_distribution<>(mu_prime, sigma));
    }
    return dists;
}

// sample from a lognormal distribution, with optional use of cutoff factor
double sample(std::lognormal_distribution<>& dist, std::mt19937& rng, bool cutoff = false) {
    double x = dist(rng);
    if (cutoff) {
        int max_tries = 100;
        const double mean = std::exp(dist.m() + 0.5 * dist.s() * dist.s());
        while (x > dist_cutoff * mean) {
            x = dist(rng);
            if (max_tries-- == 0) {
                x = mean;
                break;
            }
        }
    }
    return x;
}

// return a vector of samples from a vector of lognormal distributions
std::vector<double> sample(std::vector<std::lognormal_distribution<>>& dists, std::mt19937& rng, bool cutoff = false) {
    std::vector<double> samples;
    for (unsigned i = 0; i < dists.size(); i++) {
        samples.push_back(sample(dists[i], rng, cutoff));
    }
    return samples;
}

// print initiation message
void welcome() {
    std::cout << "--------------------------" << std::endl
              << "|  Welcome to PolyMorph  |" << std::endl
              << "--------------------------" << std::endl;
    std::cout << "Number of OpenMP threads = " << omp_get_max_threads() << std::endl;
    std::cout << "Simulation started at " << getCurrentDateTime() << std::endl;
}

// print termination message
void goodbye() {
    std::cout << "Simulation finished at " << getCurrentDateTime() << std::endl;
}

// save simluations parameters
void write_config() {
    std::ofstream config(output_folder + "/simulation.cfg");
    if (!config.is_open()) {
        std::cerr << "Error: Could not open file " << output_folder << "/simulation.cfg"<< std::endl;
        return;
    }
    config
        << "Date=" << __DATE__ << std::endl
        << "Time=" << __TIME__ << std::endl
        << "h=" << h << std::endl
        << "lmin=" << lmin << std::endl
        << "lmax=" << lmax << std::endl
        << "C=" << C << std::endl
        << "alpha_mu=" << alpha_mu << std::endl
        << "alpha_CV=" << alpha_CV << std::endl
        << "beta=" << beta << std::endl
        << "Amin=" << Amin << std::endl
        << "Amax_mu=" << Amax_mu << std::endl
        << "Amax_CV=" << Amax_CV << std::endl
        << "gam=" << gam << std::endl
        << "ka=" << ka << std::endl
        << "kl=" << kl << std::endl
        << "kb=" << kb << std::endl
        << "kr=" << kr << std::endl
        << "kh=" << kh << std::endl
        << "sh=" << sh << std::endl
        << "ss=" << ss << std::endl
        << "mu=" << mu << std::endl
        << "rho=" << rho << std::endl
        << "g=" << g << std::endl
        << "gl=" << gl << std::endl
        << "cv=" << cv << std::endl
        << "cd=" << cd << std::endl
        << "cc=" << cc << std::endl
        << "dt=" << dt << std::endl
        << "Nf=" << Nf << std::endl
        << "Ns=" << Ns << std::endl
        << "Nr=" << Nr << std::endl
        << "NUM_SPECIES=" << NUM_SPECIES << std::endl
        << "NUM_KINETIC=" << NUM_KINETIC << std::endl
        << "ADVECTION_DILUTION=" << ADVECTION_DILUTION << std::endl
        << "RNG_SEED=" << RNG_SEED << std::endl
        << "k0=" << to_string(k0) << std::endl
        << "k_mu=" << to_string(k_mu) << std::endl
        << "k_CV=" << to_string(k_CV) << std::endl
        << "D0=" << to_string(D0) << std::endl
        << "D_mu=" << to_string(D_mu) << std::endl
        << "D_CV=" << to_string(D_CV) << std::endl
        << "anisotropy=" << to_string(anisotropy) << std::endl
        << "domain=[" << domain.x0 << ", " << domain.y0 << ", " << domain.x1 << ", " << domain.y1 << "]" << std::endl
        << "dx=" << dx << std::endl
        << "dist_cutoff=" << dist_cutoff << std::endl
        << "kd=" << kd << std::endl
        << "IDW_rcut=" << IDW_rcut << std::endl
        << "input_file=" << input_file << std::endl
        << "output_folder=" << output_folder << std::endl
        << "Output::c=" << Output::c << std::endl
        << "Output::grad_c=" << Output::grad_c << std::endl
        << "Output::D=" << Output::D << std::endl
        << "Output::k=" << Output::k << std::endl
        << "Output::parent_idx=" << Output::parent_idx << std::endl
        << "Output::cell_type=" << Output::cell_type << std::endl
        << "Output::v=" << Output::v << std::endl;
    config.close();
}

#endif
