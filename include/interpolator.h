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
#include "utils.h"

enum class InterpolationMethod {
  IDW, // inverse distance weighting: good but expensive (adjust cutoff_radius in param.h)
  BILINEAR, // bilinear interpolation: cheaper but some artifacts
  ZERO // set velocity to zero outside of tissue: most efficient when only tissue interior is needed
};

// This class takes care of the data scattering and gathering between ensemble and solver
struct Interpolator {
  Ensemble& ensemble;
  Solver& solver;
  int istart, jstart, iend, jend; // limits of the grid points to be updated (within ensemble bounds)
  Grid<int>& prev_idx; // stores index of polygon in which a grid point lies. Negative indices indicate a background node.
  Grid<int> new_idx;  // temporary grid to store the new parent_idx
  InterpolationMethod ext_interpolation_method; // velocity interpolation method for exterior gridpoints (not inside polygons)
  
  Interpolator(Ensemble& ensemble, Solver& solver) : ensemble(ensemble), solver(solver), prev_idx(solver.parent_idx) {
    new_idx = Grid<int>(solver.Nx, solver.Ny, -1);
    ext_interpolation_method = InterpolationMethod::ZERO; // default method
  }
  
  // scatter coefficients D, k, p and velocity from polygons to grid points
  void scatter();

  // gather concentration c from grid points to polygons
  // important: depends on scatter being called every iteration to build the parent_idx
  void gather();

  // Search algorithm to find parent polygon for a grid point.
  std::size_t find_parent(Point grid_point);

  // returns the interpolated velocity at a grid point inside a parent polygon using IDW
  Point interior_IDW_vel_interpolation(const Point& grid_point, const int parent_idx);

  // IDW interpolation of velocity field at background nodes
  Point ext_IDW_vel_interpolation(int i, int j, double cutoff_radius);

  // bilinear interpolation of velocity field at background nodes
  Point ext_bilinear_vel_interpolation(int i, int j);
};

void Interpolator::scatter() {
  istart = std::max(int((ensemble.x0 - solver.domain.x0) / solver.dx) + 1, 0);
  jstart = std::max(int((ensemble.y0 - solver.domain.y0) / solver.dx) + 1, 0);
  iend = std::min(int((ensemble.x1 - solver.domain.x0) / solver.dx), solver.Nx);
  jend = std::min(int((ensemble.y1 - solver.domain.y0) / solver.dx), solver.Ny);

  #pragma omp parallel
  {
    #pragma omp for collapse(2)
    for (int i = istart; i < iend; i++) {
      for (int j = jstart; j < jend; j++) { 
        const double x = solver.domain.x0 + i * solver.dx;
        const double y = solver.domain.y0 + j * solver.dx;
        const Point grid_point(x, y);
        // check if still the same parent (ignore rigid)
        if (prev_idx(i, j) >= Nr && prev_idx(i, j) < (int) ensemble.polygons.size() && ensemble.polygons[prev_idx(i, j)].contains(grid_point)) {
          new_idx(i, j) = prev_idx(i, j);
        } 
        else {
          new_idx(i, j) = find_parent(grid_point);
        }
        // scatter values to grid
        if (new_idx(i, j) < Nr) { // is background node
          solver.D(i, j) = D0; // background diffusion
          solver.k(i, j) = k0; // background kinetic coefficients
        } else { 
          solver.D(i, j) = ensemble.polygons[new_idx(i, j)].D;
          solver.k(i, j) = ensemble.polygons[new_idx(i, j)].k;
          if (ADVECTION_DILUTION_EN) {
            if (0 <= i && i < solver.Nx - 1 && 0 <= j && j < solver.Ny - 1) { // leave boundary nodes at zero
              solver.velocity(i, j) = interior_IDW_vel_interpolation(grid_point, new_idx(i, j));
            }
          }
        }
      }
    }
    // update parent_idx
    solver.parent_idx.parallel_copy_from(new_idx);

    // interpolate remaining velocity field at background nodes
    if (ADVECTION_DILUTION_EN) {
      #pragma omp for collapse(2)
      for (int i = 1; i < solver.Nx - 1; i++) {
        for (int j = 1; j < solver.Ny - 1; j++) {
          if (solver.parent_idx(i, j) < 0) { // only treat real background nodes. No velocity in rigid polygons
            if (ext_interpolation_method == InterpolationMethod::BILINEAR) {
              solver.velocity(i, j) = ext_bilinear_vel_interpolation(i, j);
            } else if (ext_interpolation_method == InterpolationMethod::IDW) {
              solver.velocity(i, j) = ext_IDW_vel_interpolation(i, j, IDW_cutoff_radius);
            } else if (ext_interpolation_method == InterpolationMethod::ZERO) {
              solver.velocity(i, j) = Point(0, 0);
            } else {
              std::cerr << "Unknown interpolation method" << std::endl;
            }
          } 
        }
      }
    }
  }
}

void Interpolator::gather() {
  // get all children from parent idx built during scatter()
  for (auto& polygon : ensemble.polygons) {
    polygon.children.clear();
  }
  // cannot easily parallelize this part
  for (int i = istart; i < iend; i++) {
    for (int j = jstart; j < jend; j++) {
      if (solver.parent_idx(i, j) >= int(Nr)) { // skip background nodes
        ensemble.polygons[solver.parent_idx(i, j)].children.push_back(Index(i, j)); 
      }
    }
  }
  // accumulate average concentration and concentration gradient from children
  #pragma omp parallel for
  for (size_t p = Nr; p < ensemble.polygons.size(); p++) {
    auto& polygon = ensemble.polygons[p];
    const int num_children = polygon.children.size();
    if (num_children == 0) continue; // skip if no interior grid points
    const double inv_size = 1.0 / num_children;

    for (int sp = 0; sp < NUM_SPECIES; sp++){
      polygon.c[sp] = 0.0; // reset values
      polygon.grad_c[sp] = Point(0, 0); // reset values
      for (const Index& idx : polygon.children) {
        polygon.c[sp] += inv_size * solver.c(idx)[sp];
        polygon.grad_c[sp].add(inv_size, solver.grad_c(idx)[sp]);
      }
    }
    // Determine cell_type
    polygon.cell_type = ensemble.cellType(polygon, polygon.c, polygon.grad_c, ensemble.t);
  } 
}

std::size_t Interpolator::find_parent(Point grid_point) {
  std::size_t bxi = (grid_point.x - ensemble.x0) / ensemble.bs + 1; // box index in x direction
  const std::size_t byi = (grid_point.y - ensemble.y0) / ensemble.bs + 1; // box index in y direction
  std::unordered_set<std::size_t> checked_polygons; // store checked polygons to avoid checking them again
  bool last_iteration = false; // abort search as soon as we can be sure it's a background node
  while (bxi < ensemble.Nx) {
    for (Vertex* v = ensemble.first[bxi * ensemble.Ny + byi]; v; v = v->next) { // loop over vertices in box
      if (checked_polygons.find(v->p) == checked_polygons.end()) {  // new polygon encountered
        if (v->p >= Nr && ensemble.polygons[v->p].contains(grid_point)) { // don't want to check rigid polygons
          return v->p;
        } else {
          checked_polygons.insert(v->p); // maybe more efficient to omit the set. not many vert of same polygon in box
        }
      }
    }
    if (last_iteration) return -2; // background node
    if (!checked_polygons.empty()) last_iteration = true; // only go 1 more layer
    ++bxi;
  }
  return -2; // background node (reached boundary of ensemble box)
}

Point Interpolator::interior_IDW_vel_interpolation(const Point& grid_point, const int parent_idx) {
  const Polygon& parent = ensemble.polygons[parent_idx];
  std::vector<double> weights;
  double total_weight = 0;
  for (size_t i = 0; i < parent.vertices.size(); i++) {
    const Vertex& vertex = parent.vertices[i];
    double distance = (vertex.r - grid_point).length();
    double weight = 1.0 / (distance + 1e-6 * h); // avoid division by zero
    weights.push_back(weight);
    total_weight += weight;
  }    
  Point vel = Point(0, 0);
  for (size_t i = 0; i < parent.vertices.size(); i++) {
    const Vertex& vertex = parent.vertices[i];
    vel = vel + weights[i] / total_weight * vertex.v;
  }
  return vel;
}

Point Interpolator::ext_bilinear_vel_interpolation(int i, int j) {
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
  const double eps = 1e-6 * h; // avoid division by zero
  const double w_left = dx_right / (dx_left + dx_right + eps);
  const double w_right = dx_left / (dx_left + dx_right + eps);
  const double w_down = dy_up / (dy_down + dy_up + eps);
  const double w_up = dy_down / (dy_down + dy_up + eps);
  // interpolate velocity
  Point vel;
  vel.x = w_left * solver.velocity(i_left, j).x + w_right * solver.velocity(i_right, j).x; 
  vel.y = w_down * solver.velocity(i, j_down).y + w_up * solver.velocity(i, j_up).y;
  return vel;
}

Point Interpolator::ext_IDW_vel_interpolation(int i, int j, double cutoff_radius) { // maybe interpolate with vertices inside the local box
  int cutoff_index = cutoff_radius / solver.dx;
  double total_weight = 1e-6 * h; // avoid division by zero
  Point velocity = Point(0, 0);
  for (int ii = i - cutoff_index; ii <= i + cutoff_index; ii++) {
    for (int jj = j - cutoff_index; jj <= j + cutoff_index; jj++) {
      double distance2 = (i - ii) * (i - ii) + (j - jj) * (j - jj);
      // skip nodes outside cutoff radius
      if (distance2 > cutoff_index * cutoff_index) continue; 
      // ensure indices are within bounds
      if (ii >= 0 && ii < solver.Nx && jj >= 0 && jj < solver.Ny) {
        // use nodes inside polygon and boundary nodes
        if (solver.parent_idx(ii, jj) >= 0 || ii == 0 || ii == solver.Nx-1 || jj == 0 || jj == solver.Ny-1) { 
          double weight = 1.0 / (distance2 + 1e-6 * h); // avoid division by zero
          velocity.add(weight, solver.velocity(ii, jj)); 
          total_weight += weight;
        }
      }
    }
  }
  return 1.0 / total_weight * velocity;
}

#endif