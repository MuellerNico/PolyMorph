// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <vector>

// 2D vector
struct Point {
    double x, y; // coordinates
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
    // lexicographical comparison for std::map (used in Ensemble::write_OFF)
    bool operator<(const Point& r) const {
        if (x == r.x) return y < r.y;
        return x < r.x;
    }
};
Point operator*(const double a, const Point& r) { return {a * r.x, a * r.y}; }

// grid indices in x and y direction
struct Index { int i, j; };

// rectangular spatial domain
struct Domain { double x0, y0, x1, y1; };

// shorted squared distance from point to edge
double point_edge_dist2(const Point& r0, const Point& r1, const Point& r2, double& xi, double& xit, Point& dr) {
    const Point r12 = r2 - r1;
    const Point r10 = r0 - r1;
    xi = r10 * r12 / (r12 * r12);
    xit = std::min(std::max(xi, 0.), 1.);
    dr = r10 - xit * r12;
    return dr * dr;
}

// polygon vertex
struct Vertex {
    Point r, v, a; // position, velocity, acceleration
    std::size_t p; // polygon index
    Vertex* next; // pointer to next vertex in same box
    double l0; // rest length of edge to the right
};

// polygonal cell
struct Polygon {
    std::vector<Vertex> vertices; // vertex list in counter-clockwise orientation
    bool phase; // phase of the enclosed medium
    double A0, A, Amax, alpha; // target, actual & division area, area growth rate
    std::vector<double> D, k; // diffusion coefficients, kinetic coefficients
    std::vector<double> c; // local concentrations
    std::vector<Point> grad_c; // local concentration gradients
    std::vector<Index> children; // stores the x,y indices of grid points (nodes) that lie inside the polygon
    int cell_type = 0; // cell type (or general-purpose label for different applications)

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
    // checks whether a point lies inside this polygon
    bool contains(const Point& r) const {
        bool in = false;
        for (std::size_t i = vertices.size() - 1, j = 0; j < vertices.size(); i = j++)
            if ((vertices[i].r.y > r.y) != (vertices[j].r.y > r.y) &&
                r.x < (vertices[i].r.x - vertices[j].r.x) * (r.y - vertices[j].r.y) / (vertices[i].r.y - vertices[j].r.y) + vertices[j].r.x)
                in = !in;
        return in;
    }

    // index of the polygon in the ensemble
    std::size_t index() const {
        return vertices[0].p; // vertices store the global polygon index they belong to
    }
};

#endif // GEOMETRY_H
