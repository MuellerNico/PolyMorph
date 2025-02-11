#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "utils.h"

struct Point { // basically a 2D vector
  double x, y; // coordinates
  Point () : x(0), y(0) {} // default constructor
  Point (const double x, const double y) : x(x), y(y) {} // constructor
  Point operator+(const Point& r) const { return {x + r.x, y + r.y}; } // vector sum
  Point operator-(const Point& r) const { return {x - r.x, y - r.y}; } // vector difference
  Point cross() const { return {-y, x}; } // perpendicular vector
  double operator*(const Point& r) const { return x * r.x + y * r.y; } // dot product
  double operator^(const Point& r) const { return x * r.y - y * r.x; } // wedge product
  double length2() const { return x * x + y * y; } // squared norm
  double length() const { return std::sqrt(length2()); } // norm
  // multiply and add (+= a*r)
  void add(const double a, const Point& r) {
    #pragma omp atomic
    x += a * r.x;
    #pragma omp atomic
    y += a * r.y;
  }
  // lexographical comparison for std::map (used in write_OFF)
  bool operator<(const Point& r) const {
    if (x == r.x) return y < r.y;
    return x < r.x;
  }
};

// scalar times vector multiplication
Point operator*(const double a, const Point& r) { return {a * r.x, a * r.y}; }

double point_edge_dist2(const Point& r0, const Point& r1, const Point& r2, double& xi, double& xit, Point& dr) {
  const Point r12 = r2 - r1;
  const Point r10 = r0 - r1;
  xi = r10 * r12 / (r12 * r12);
  xit = std::min(std::max(xi, 0.), 1.);
  dr = r10 - xit * r12;
  return dr * dr;
}

struct Vertex {
  Point r, v, a; // position, velocity, acceleration
  std::size_t p; // polygon index
  Vertex* next; // pointer to next vertex in same box
  double l0; // rest length of edge to the right
};

struct Polygon {
  std::vector<Vertex> vertices; // vertex list in counter-clockwise orientation
  bool phase; // phase of the enclosed medium
  double A0, A, Amax, alpha; // target, actual & division area, area growth rate
  std::vector<double> D, k, p; // diffusion constant, kinetic coefficients, production rate
  std::vector<double> c = std::vector<double>(NUM_SPECIES, 0); // local concentration
  std::vector<Point> grad_c = std::vector<Point>(NUM_SPECIES, Point(0, 0)); // local concentration gradient
  std::vector<Index> children; // stores the indices of grid points (nodes) that lie inside the polygon
  int cell_type = 0; // cell type (or general purpose flag for different user applications)

  double area() {
    A = 0;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      A += vertices[i].r ^ vertices[j].r;
    return A /= 2;
  }
  double perimeter() const {
    double L = 0;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      L += (vertices[j].r - vertices[i].r).length();
    return L;
  }
  double perimeter0() const {
    double L0 = 0;
    for (auto& v : vertices)
      L0 += v.l0;
    return L0;
  }
  // returns true if the point lies inside this polygon
  bool contains(const Point& r) const {
    bool in = false;
    for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
      if ((vertices[i].r.y > r.y) != (vertices[j].r.y > r.y) &&
          r.x < (vertices[i].r.x - vertices[j].r.x) * (r.y - vertices[j].r.y) / (vertices[i].r.y - vertices[j].r.y) + vertices[j].r.x)
        in = !in;
    return in; 
  }
  // calculate center of polygon
  Point midpoint() const {
    Point midp;
    for (const Vertex& v : vertices) {
      midp = midp + v.r;
    }
    return 1.0/vertices.size() * midp;
  }

  // returns the index of the polygon in the global ensemble list
  size_t polygon_index() const {
    return vertices[0].p; // vertices store the global polygon index they belong to
  }
};

#endif // GEOMETRY_H