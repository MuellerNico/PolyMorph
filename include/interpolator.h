// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef INTERPOLATOR_H
#define INTERPOLATOR_H

#include <unordered_set>

#include "solver.h"
#include "ensemble.h"

enum class InterpolationMethod {
    IDW, // inverse distance weighting: good but expensive (adjust cutoff_radius in param.h)
    BILINEAR, // bilinear interpolation: cheaper but some artifacts
    ZERO // set velocity to zero outside of tissue: most efficient when only tissue interior is needed
};

// class taking care of data scattering and gathering between ensemble and solver
struct Interpolator {
    Ensemble& ensemble;
    Solver& solver;
    int istart, jstart, iend, jend; // limits of the grid points to be updated (within ensemble bounds)
    Grid<int> new_idx; // temporary grid to store the new parent_idx
    InterpolationMethod exterior_method; // velocity interpolation method for exterior gridpoints (not inside polygons)
    static constexpr double eps = 1e-6 * h; // short lengthscale to void division by zero

    Interpolator(Ensemble& ensemble, Solver& solver) : ensemble(ensemble), solver(solver),
                 new_idx(solver.Nx, solver.Ny, -1), exterior_method(InterpolationMethod::ZERO) {}

    // scatter coefficients D, k and velocity from polygons to grid points
    void scatter();

    // gather concentration c from grid points to polygons
    // important: depends on scatter being called every iteration to build the parent_idx
    void gather();

    // search algorithm to find parent polygon for a grid point
    int find_parent(Point grid_point);

    // inverse distance weighting (IDW) interpolation of velocity field at a node inside a parent polygon
    Point interior_IDW(const Point& grid_point, const Polygon& p);

    // IDW interpolation of velocity field at background nodes
    Point exterior_IDW(int i, int j);

    // bilinear interpolation of velocity field at background nodes
    Point exterior_bilinear(int i, int j);
};

void Interpolator::scatter() {
    istart = std::max(static_cast<int>((ensemble.bbox.x0 - domain.x0) / dx) + 1, 0);
    jstart = std::max(static_cast<int>((ensemble.bbox.y0 - domain.y0) / dx) + 1, 0);
    iend   = std::min(static_cast<int>((ensemble.bbox.x1 - domain.x0) / dx), solver.Nx);
    jend   = std::min(static_cast<int>((ensemble.bbox.y1 - domain.y0) / dx), solver.Ny);
  
    Grid<int>& prev_idx = solver.parent_idx; // index of polygon in which a grid point lies (negative indices -1 or -2 indicate background)

    #pragma omp parallel
    {
        // fill new_idx with background
        #pragma omp for collapse(2)
        for (int i = 0; i < solver.Nx; i++) {
            for (int j = 0; j < solver.Ny; j++) {
                new_idx(i, j) = -1;
            }
        }
        
        #pragma omp for collapse(2)
        for (int i = istart; i < iend; i++) {
            for (int j = jstart; j < jend; j++) {
                const double x = domain.x0 + i * dx;
                const double y = domain.y0 + j * dx;
                const Point grid_point{x, y};
                // check if still the same parent (ignore rigid polygons)
                if (prev_idx(i, j) >= static_cast<int>(Nr) && prev_idx(i, j) < static_cast<int>(ensemble.polygons.size()) && ensemble.polygons[prev_idx(i, j)].contains(grid_point)) {
                    new_idx(i, j) = prev_idx(i, j);
                } else {
                    new_idx(i, j) = find_parent(grid_point);
                }
                // scatter values to grid
                if (new_idx(i, j) < static_cast<int>(Nr)) { // is background node or in a rigid polygon
                    solver.D(i, j) = D0; // background diffusion
                    solver.k(i, j) = k0; // background kinetic coefficients
                } else {
                    solver.D(i, j) = ensemble.polygons[new_idx(i, j)].D;
                    solver.k(i, j) = ensemble.polygons[new_idx(i, j)].k;
                    if (ADVECTION_DILUTION) {
                        if (0 <= i && i < solver.Nx - 1 && 0 <= j && j < solver.Ny - 1) { // leave boundary nodes at zero
                            solver.v(i, j) = interior_IDW(grid_point, ensemble.polygons[new_idx(i, j)]);
                        }
                    }
                }
            }
        }
        
        // interpolate remaining velocity field at background nodes
        if (ADVECTION_DILUTION) {
            #pragma omp for collapse(2)
            for (int i = 1; i < solver.Nx - 1; i++) {
                for (int j = 1; j < solver.Ny - 1; j++) {
                    if (new_idx(i, j) < 0) { // only treat real background nodes (no velocity in rigid polygons)
                        if (exterior_method == InterpolationMethod::BILINEAR) {
                            solver.v(i, j) = exterior_bilinear(i, j);
                        } else if (exterior_method == InterpolationMethod::IDW) {
                            solver.v(i, j) = exterior_IDW(i, j);
                        } else if (exterior_method == InterpolationMethod::ZERO) {
                            solver.v(i, j) = {0, 0};
                        } else {
                            std::cerr << "Error: Unknown interpolation method" << std::endl;
                        }
                    }
                }
            }
        }
    }

    prev_idx.swap(new_idx);
}

void Interpolator::gather() {
    // get all children from parent index built during scatter()
    for (auto& polygon : ensemble.polygons) {
        polygon.children.clear();
    }
    // cannot easily be parallelized
    for (int i = istart; i < iend; i++) {
        for (int j = jstart; j < jend; j++) {
          if (solver.parent_idx(i, j) >= static_cast<int>(Nr)) { // skip background nodes
            ensemble.polygons[solver.parent_idx(i, j)].children.push_back({i, j});
            }
        }
    }
    // accumulate average concentration and concentration gradient from children
    #pragma omp parallel for
    for (std::size_t p = Nr; p < ensemble.polygons.size(); p++) {
        Polygon& polygon = ensemble.polygons[p];
        if (polygon.children.size() == 0) continue; // skip if no interior grid points
        const double inv_size = 1.0 / polygon.children.size();
        for (std::size_t sp = 0; sp < NUM_SPECIES; sp++) {
            polygon.c[sp] = 0; // reset values
            polygon.grad_c[sp] = {0, 0}; // reset values
            for (const Index& idx : polygon.children) {
                polygon.c[sp] += inv_size * solver.c(idx.i, idx.j)[sp];
                polygon.grad_c[sp].add(inv_size, solver.grad_c(idx.i, idx.j)[sp]);
            }
        }
        polygon.cell_type = ensemble.cellType(polygon, polygon.c, polygon.grad_c, ensemble.t);
    }
}

int Interpolator::find_parent(Point grid_point) {
          std::size_t bxi = (grid_point.x - ensemble.bbox.x0) / ensemble.bs + 1; // box index in x direction
    const std::size_t byi = (grid_point.y - ensemble.bbox.y0) / ensemble.bs + 1; // box index in y direction
    std::unordered_set<std::size_t> checked_polygons; // store checked polygons to avoid checking them again
    bool last_iteration = false; // abort search as soon as it's a background node for sure
    for (; bxi < ensemble.Nx; bxi++) {
        for (Vertex* v = ensemble.first[bxi * ensemble.Ny + byi]; v; v = v->next) { // loop over vertices in box
            if (checked_polygons.find(v->p) == checked_polygons.end()) { // new polygon encountered
                if (v->p >= Nr && ensemble.polygons[v->p].contains(grid_point)) { // skip rigid polygons
                    return v->p;
                } else {
                    checked_polygons.insert(v->p); // may be more efficient to omit the set (not many vertices of same polygon in box)
                }
            }
        }
        if (last_iteration) return -2; // background node
        if (!checked_polygons.empty()) last_iteration = true; // only go 1 more layer
    }
    return -2; // background node (reached boundary of ensemble box)
}

Point Interpolator::interior_IDW(const Point& grid_point, const Polygon& p) {
    std::vector<double> weights(p.vertices.size());
    double total_weight = 0;
    for (std::size_t i = 0; i < p.vertices.size(); i++) {
        const double distance = (p.vertices[i].r - grid_point).length();
        const double weight = 1 / (distance + eps);
        weights[i] = weight;
        total_weight += weight;
    }
    Point vel{0, 0};
    for (std::size_t i = 0; i < p.vertices.size(); i++) {
        vel.add(weights[i] / total_weight, p.vertices[i].v);
    }
    return vel;
}

Point Interpolator::exterior_bilinear(int i, int j) {
    // find first non-background node or boundary node in each direction
    int i_left = i - 1;
    while (i_left > 0 && solver.parent_idx(i_left, j) < 0) {
        i_left--;
    }
    int i_right = i + 1;
    while (i_right < solver.Nx-1 && solver.parent_idx(i_right, j) < 0) {
        i_right++;
    }
    int j_down = j - 1;
    while (j_down > 0 && solver.parent_idx(i, j_down) < 0) {
        j_down--;
    }
    int j_up = j + 1;
    while (j_up < solver.Ny-1 && solver.parent_idx(i, j_up) < 0) {
        j_up++;
    }
    // calculate distances
    const double dx_left = i - i_left;
    const double dx_right = i_right - i;
    const double dy_down = j - j_down;
    const double dy_up = j_up - j;
    // calculate weights
    const double w_left = dx_right / (dx_left + dx_right + eps);
    const double w_right = dx_left / (dx_left + dx_right + eps);
    const double w_down = dy_up / (dy_down + dy_up + eps);
    const double w_up = dy_down / (dy_down + dy_up + eps);
    // interpolate velocity
    Point vel;
    vel.x = w_left * solver.v(i_left, j).x + w_right * solver.v(i_right, j).x;
    vel.y = w_down * solver.v(i, j_down).y + w_up    * solver.v(i,    j_up).y;
    return vel;
}

Point Interpolator::exterior_IDW(int i, int j) {
    static const int cutoff_index = IDW_rcut / dx;
    const int iistart = std::max(i - cutoff_index, 0);
    const int jjstart = std::max(j - cutoff_index, 0);
    const int iiend   = std::min(i + cutoff_index, solver.Nx - 1);
    const int jjend   = std::min(j + cutoff_index, solver.Ny - 1);
    double total_weight = eps; // avoid division by zero
    Point vel{0, 0};
    for (int ii = iistart; ii <= iiend; ii++) {
        for (int jj = jjstart; jj <= jjend; jj++) {
            const double distance2 = (i - ii) * (i - ii) + (j - jj) * (j - jj);
            // skip nodes outside cutoff radius
            if (distance2 <= cutoff_index * cutoff_index) {
                // use nodes inside polygon and boundary nodes
                if (solver.parent_idx(ii, jj) >= 0 || ii == 0 || ii == solver.Nx-1 || jj == 0 || jj == solver.Ny-1) {
                    const double weight = 1 / (distance2 + eps * eps);
                    vel.add(weight, solver.v(ii, jj));
                    total_weight += weight;
                }
            }
        }
    }
    return 1 / total_weight * vel;
}

#endif
