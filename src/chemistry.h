#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include <numeric>

#include "polyhoop.h"
#include "utils.h"

// Handles polygon modification due to signaling effects. 
// Could also be done in ensemble.step()
struct Chemistry {
  Ensemble& ensemble;
  bool growth_control = false; // set true to stop growth when below threshold
  std::function<bool(Polygon&)> is_producing = [](Polygon& p) { return false;};

  Chemistry(Ensemble& ensemble) : ensemble(ensemble) {}

  void update() { // ToDo: split into flag() and produce()
    #pragma omp parallel for
    for (int i = Nr; i < ensemble.polygons.size(); i++) {
      auto& cell = ensemble.polygons[i];
      // set production
      double p_sum = std::accumulate(cell.p.begin(), cell.p.end(), 0.0); // assuming p_i >= 0. ToDo: allow sinks
      if (p_sum == 0 && is_producing(cell)) {
        cell.p = sample(p_dist, rng);
      } else if (p_sum > 0 && !is_producing(cell)) {
        cell.p = std::vector<double>(NUM_SPECIES, 0);
      }
      // flag below threshold
      if (!cell.flag && cell.u[0] < cell.threshold[0]) {  // atm only consider first species
        cell.flag = true;
        if (growth_control) cell.alpha = 0;
      } else if (cell.flag && cell.u[0] > cell.threshold[0]){
        cell.flag = false;
        if (growth_control) cell.alpha = cell.alpha0; 
      }
    }
  }  

  double get_border_sharpness() { // "width" of border
    double xmin = ensemble.x1;
    double xmax = ensemble.x0;
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      for (auto& vertex : cell.vertices) {
        if (cell.flag) {
          xmin = std::min(xmin, vertex.r.x); // leftmost flagged point
        } else {
          xmax = std::max(xmax, vertex.r.x); // rightmost unflagged
        }
      }
    }
    std::cout << "xmin=" << xmin << " xmax=" << xmax << std::endl;
    return xmax - xmin;
  }

  void chemotaxis() {
    #pragma omp parallel for
    for (int p = Nr; p < ensemble.polygons.size(); p++) {
      auto& cell = ensemble.polygons[p];
      for (auto& vertex : cell.vertices) {
        // move towards higher concentration
        // manipulate vertex acceleration
        // Q: how to get morphogen gradient cheaply?
        ; 
      }
    }
  }
};

#endif