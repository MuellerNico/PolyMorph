// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <random>
#include <map>
#include <cmath>
#include <functional>

#include "param.h"
#include "geometry.h"
#include "utils.h"

// User-defined lambda functions to control concentration effects on cell behavior
template <typename T>
using Function = std::function<T(const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t)>;

// User-defined lambda functions to control concentration effects on cell contact
template <typename T>
using PairFunction = std::function<T(const Polygon& p1, const Polygon& p2, const std::vector<double>& c1, const std::vector<double>& c2, const std::vector<Point>& grad_c1, const std::vector<Point>& grad_c2, double t)>;

struct Ensemble {
    std::vector<Polygon> polygons; // global list of polygons
    std::vector<Vertex*> first; // pointer to first vertex in each box
    std::size_t Nx, Ny; // number of boxes in x,y direction
    Domain bbox; // bounding box
    double bs; // box size
    double t; // time
    
    std::mt19937 rng; // random number generator

    // probability distributions
    const double Amax_lnCV = std::log(1 + Amax_CV*Amax_CV);
    const double alpha_lnCV = std::log(1 + alpha_CV*alpha_CV);
    std::lognormal_distribution<> Amax_dist = std::lognormal_distribution<>(std::log(Amax_mu) - Amax_lnCV/2, std::sqrt(Amax_lnCV)); // division area distribution
    std::lognormal_distribution<> alpha_dist = std::lognormal_distribution<>(std::log(alpha_mu) - alpha_lnCV/2, std::sqrt(alpha_lnCV)); // area growth rate distribution
    std::uniform_real_distribution<> uni_dist;

    std::vector<std::lognormal_distribution<>> D_dist = create_lognormal(D_mu, D_CV);
    std::vector<std::lognormal_distribution<>> k_dist = create_lognormal(k_mu, k_CV);

    // lambda functions
    // concentration effect on single-cell properties
    Function<Point> acceleration = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return Point{0, 0}; };
    Function<int> cellType = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return 0; };
    Function<double> growthRate = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return p.alpha; };
    Function<double> maxArea = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return p.Amax; };
    Function<double> areaStiffness = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return ka; };
    Function<double> lineTension = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return gam; };
    Function<double> edgeContractilityStiffness = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return kl; };
    Function<double> bendingStiffness = [](const Polygon& p, const std::vector<double>& c, const std::vector<Point>& grad_c, double t) { return kb; };
    // concentration effect on cell-cell contact properties
    PairFunction<double> adhesionStiffness = [](const Polygon& p1, const Polygon& p2, const std::vector<double>& c1, const std::vector<double>& c2, const std::vector<Point>& grad_c1, const std::vector<Point>& grad_c2, double t) { return kh; };
    PairFunction<double> dynamicFriction = [](const Polygon& p1, const Polygon& p2, const std::vector<double>& c1, const std::vector<double>& c2, const std::vector<Point>& grad_c1, const std::vector<Point>& grad_c2, double t) { return mu; };

    Ensemble(const char* name, int seed = RNG_SEED) : t(0) {
        // initialize random number generator
        rng.seed(seed);
        
        // read OFF file header
        std::ifstream file(name);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << name << std::endl;
            std::exit(1);
        }
        file.ignore(3); // ignore "OFF"
        std::size_t Nv, Np, Ne; // number of vertices, polygons, edges
        file >> Nv >> Np >> Ne; // Ne unused

        // read all vertices
        std::vector<Point> points(Nv);
        std::vector<int> z(Nv);
        for (std::size_t i = 0; i < Nv; ++i)
            file >> points[i].x >> points[i].y >> z[i];
        
        // read all polygons
        polygons.resize(Np);
        for (std::size_t p = 0, j; p < Np; ++p) {
            file >> Nv;
            for (std::size_t i = 0; i < Nv; ++i) {
                file >> j;
                polygons[p].vertices.push_back({points[j], {0, 0}, {0, 0}, p});
            }
            polygons[p].phase = std::abs(z[j]) % 2; // the z coordinate is used as the phase
            polygons[p].A0 = polygons[p].area(); // use the initial area as target area
            polygons[p].Amax = sample(Amax_dist, rng);
            polygons[p].alpha = sample(alpha_dist, rng);

            polygons[p].D = sample(D_dist, rng, true);
            polygons[p].k = sample(k_dist, rng);
            polygons[p].c = std::vector<double>(NUM_SPECIES, 0);
            polygons[p].grad_c = std::vector<Point>(NUM_SPECIES, {0, 0});
          
            for (std::size_t i = Nv - 1, j = 0; j < Nv; i = j++)
                polygons[p].vertices[i].l0 = (polygons[p].vertices[j].r - polygons[p].vertices[i].r).length(); // edge rest length
        }
    }
    
    void remove(std::size_t p) {
        polygons[p] = polygons.back(); // copy the last polygon to index p
        for (auto& v : polygons[p].vertices)
            v.p = p; // let the vertices know about their new position
        polygons.pop_back(); // remove the last polygon
    }
    
    void step(double dt) {
        // polygon removal and division
        for (std::size_t p = Nr; p < polygons.size(); ++p) {
            auto& v = polygons[p].vertices;
            double Amax_new = maxArea(polygons[p], polygons[p].c, polygons[p].grad_c, t);

            if (polygons[p].A < Amin || v.size() < 3) {
                remove(p--);
            }
            else if (polygons[p].A > Amax_new && v.size() > 5) {
                // compute polygon centroid, inertia tensor and division axis
                Point c{0, 0};
                double Ixx = 0, Iyy = 0, Ixy = 0;
                for (std::size_t i = v.size() - 1, j = 0; j < v.size(); i = j++) {
                    const Point rsum = v[i].r + v[j].r;
                    const double w = v[i].r ^ v[j].r;
                    c.add(w, rsum);
                    Ixx += w * (rsum.y * rsum.y - v[i].r.y * v[j].r.y);
                    Iyy += w * (rsum.x * rsum.x - v[i].r.x * v[j].r.x);
                    Ixy += w * (rsum.x * rsum.y + v[i].r.x * v[i].r.y + v[j].r.x * v[j].r.y);
                }
                c = 1 / (6 * polygons[p].A) * c;
                Ixx =  Ixx / 12 - c.y * c.y * polygons[p].A;
                Iyy =  Iyy / 12 - c.x * c.x * polygons[p].A;
                Ixy = -Ixy / 24 + c.x * c.y * polygons[p].A;
                const double dI = (Ixx - Iyy) / 2;
                const double eig = dI + Iyy + std::sqrt(dI * dI + Ixy * Ixy);
                Point axis = (Ixx < Iyy ? Point{Ixy, eig - Ixx} : Point{eig - Iyy, Ixy});
                
                // use a random cell division axis if the inertia tensor is proportional to the identity matrix
                const double l = axis.length();
                if (l == 0) {
                    const double angle = 2 * M_PI * uni_dist(rng);
                    axis = {std::cos(angle), std::sin(angle)};
                }
                else
                    axis = 1 / l * axis; // otherwise, normalize the axis
                
                // divide the polygon in half, ensuring the end vertices are separated enough
                std::size_t vend [4];
                for (unsigned int a = 0, i = v.size() - 1, j = 0; a < 4; i = j++) {
                    const double d0 = axis ^ (v[i].r - c);
                    const double d1 = axis ^ (v[j].r - c);
                    if (d0 * d1 <= 0) {
                        vend[a] = i;
                        vend[a+1] = j;
                        for (int d = 0; (v[vend[a+1]].r - v[vend[a]].r).length2() < h * h; d = 1 - d)
                            vend[a+d] = (vend[a+d] + v.size() + 2*d-1) % v.size();
                        j = (vend[a+1] + 1) % v.size();
                        a += 2;
                    }
                }
                
                // make sure the new edges are stretched just as much as the rest of the polygon
                const double ls = polygons[p].perimeter() / polygons[p].perimeter0(); // average stretch ratio
                v[vend[0]].l0 = (v[vend[3]].r - v[vend[0]].r).length() / ls;
                v[vend[2]].l0 = (v[vend[1]].r - v[vend[2]].r).length() / ls;
                
                // distribute the vertices to the two new polygons
                const double A0 = polygons[p].A0;
                std::vector<Vertex> vold;
                vold.swap(v);
                polygons[p] = {{vold[vend[2]]}, polygons[p].phase, 0, 0, polygons[p].Amax, polygons[p].alpha, polygons[p].D, polygons[p].k}; // new polygon 1 (inherit properties)
                polygons.push_back({{vold[vend[0]]}, polygons[p].phase, 0, 0, sample(Amax_dist, rng), sample(alpha_dist, rng), sample(D_dist, rng, true), sample(k_dist, rng)}); // new polygon 2 (sample new properties)
                for (std::size_t i = vend[1]; i != vend[2]; i = (i + 1) % vold.size())
                    polygons[p].vertices.push_back(vold[i]);
                for (std::size_t i = vend[3]; i != vend[0]; i = (i + 1) % vold.size())
                    polygons.back().vertices.push_back(vold[i]);
                for (auto& v2 : polygons.back().vertices)
                    v2.p = polygons.size() - 1;
                
                // update polygon areas, splitting the target area in proportion to the actual area
                polygons[p].area();
                polygons.back().area();
                polygons[p].A0 = A0 * polygons[p].A / (polygons[p].A + polygons.back().A);
                polygons.back().A0 = A0 - polygons[p].A0;
                
                // allocate memory for the concentrations
                polygons[p].c = std::vector<double>(NUM_SPECIES, 0);
                polygons[p].grad_c = std::vector<Point>(NUM_SPECIES, {0, 0});
                polygons.back().c = std::vector<double>(NUM_SPECIES, 0);
                polygons.back().grad_c = std::vector<Point>(NUM_SPECIES, {0, 0});
            }
        }
        
        // handle vertices (refine, coarsen, compute accelerations)
        #pragma omp parallel for
        for (std::size_t p = 0; p < polygons.size(); ++p) {
            auto& v = polygons[p].vertices;
            
            // refine or coarsen the polygon
            for (std::size_t i = v.size() - 1, k = i - 1, j = 0; j < v.size(); ) {
                const double l2 = (v[j].r - v[i].r).length2();
                if (l2 > lmax * lmax) { // refine long edges, including rigid polygons
                    // insert a new vertex in the middle
                    const Point rnew = 0.5 * (v[i].r + v[j].r);
                    const Point vnew = 0.5 * (v[i].v + v[j].v);
                    v[i].l0 /= 2;
                    v.insert(v.begin() + j, {rnew, vnew, {0, 0}, p, 0, v[i].l0});
                    if (k > j) ++k;
                    if (i > j) ++i;
                }
                else if (l2 < lmin * lmin && p >= Nr && v.size() > 3) { // coarsen only non-rigid polygons with more than 3 vertices
                    // merge the two vertices
                    v[k].l0 += v[i].l0 / 2;
                    v[i].l0 = v[i].l0 / 2 + v[j].l0;
                    v[i].r = 0.5 * (v[i].r + v[j].r);
                    v[i].v = 0.5 * (v[i].v + v[j].v);
                    v.erase(v.begin() + j);
                    if (k > j) --k;
                    if (i > j) --i;
                }
                else // continue checking the next edge
                    k = i, i = j++;
            }
            
            // compute vertex accelerations
            polygons[p].area(); // update the polygon area
            const double ls = polygons[p].perimeter0() / std::sqrt(4 * M_PI * polygons[p].A0 / C); // inverse stretch ratio
            const Point a = acceleration(polygons[p], polygons[p].c, polygons[p].grad_c, t);
            for (std::size_t i = v.size() - 1, k = i - 1, j = 0; j < v.size(); k = i, i = j++) {
                const Point e1 = v[i].r - v[k].r;
                const Point e2 = v[j].r - v[i].r;
                const double l1 = e1.length();
                const double l2 = e2.length();
                const Point n = (e1 + e2).cross(); // unnormalized inward normal vector
             
                const double ka_new = areaStiffness(polygons[p], polygons[p].c, polygons[p].grad_c, t);
                const double gam_new = lineTension(polygons[p], polygons[p].c, polygons[p].grad_c, t);
                const double kl_new = edgeContractilityStiffness(polygons[p], polygons[p].c, polygons[p].grad_c, t);
                const double kb_new = bendingStiffness(polygons[p], polygons[p].c, polygons[p].grad_c, t);

                v[i].a = {0, -gl}; // edge gravitational acceleration
                v[i].a.add(ka_new / 2 * (polygons[p].A - polygons[p].A0), n); // area compressibility
                v[i].a.add(-kl_new, (ls / v[k].l0 - 1 / l1) * e1 - (ls / v[i].l0 - 1 / l2) * e2); // edge contractility
                v[i].a.add(-gam_new, (1 / l1) * e1 - (1 / l2) * e2); // line tension
                v[i].a.add(rho * g / 6 * (2*polygons[p].phase-1), (v[k].r.y + v[i].r.y + v[j].r.y) * n - Point{0, e1 ^ e2}); // hydrostatic pressure
                v[i].a.add(-cv - rho * cd / 4 * std::abs(v[i].v * n), v[i].v); // viscous damping and drag

                // bending
                const Point e0 = v[k].r - v[(k + v.size() - 1) % v.size()].r;
                const Point e3 = v[(j + 1) % v.size()].r - v[j].r;
                const double l0 = e0.length();
                const double l3 = e3.length();
                const double b0 = l0 * l1 + e0 * e1;
                const double b1 = l1 * l2 + e1 * e2;
                const double b2 = l2 * l3 + e2 * e3;
                const double a0 = (e0 ^ e1) / b0;
                const double a1 = (e1 ^ e2) / b1;
                const double a2 = (e2 ^ e3) / b2;
                v[i].a.add(-8 * kb_new, (a0 / ((l0 + l1) * b0)) * (e0.cross() - a0 * e0)
                                      - (a1 / ((l1 + l2) * b1)) * (n - a1 * (e1 - e2))
                                      + (a2 / ((l2 + l3) * b2)) * (e3.cross() + a2 * e3));

                // concentration effect on cell acceleration (same for all vertices of a polygon)
                v[i].a.add(1, a);

                // domain boundaries
                v[i].a.add((polygons[p].vertices[i].r.x < domain.x0) * kd, {(domain.x0 - polygons[p].vertices[i].r.x), 0});
                v[i].a.add((polygons[p].vertices[i].r.x > domain.x1) * kd, {(domain.x1 - polygons[p].vertices[i].r.x), 0});
                v[i].a.add((polygons[p].vertices[i].r.y < domain.y0) * kd, {0, (domain.y0 - polygons[p].vertices[i].r.y)});
                v[i].a.add((polygons[p].vertices[i].r.y > domain.y1) * kd, {0, (domain.y1 - polygons[p].vertices[i].r.y)});
            }
        }
        
        // polygon-polygon interaction
        boxes(); // place all vertices into boxes
        #pragma omp parallel for
        for (std::size_t p = Nr; p < polygons.size(); ++p) {
            for (std::size_t i = polygons[p].vertices.size() - 1, k = i - 1, j = 0; j < polygons[p].vertices.size(); k = i, i = j++) {
                const std::size_t bxi = (polygons[p].vertices[i].r.x - bbox.x0) / bs + 1; // box index in x direction
                const std::size_t byi = (polygons[p].vertices[i].r.y - bbox.y0) / bs + 1; // box index in y direction
                for (std::size_t bxj = bxi - 1; bxj <= bxi + 1; ++bxj)    // interaction within 1 box in each direction
                    for (std::size_t byj = byi - 1; byj <= byi + 1; ++byj)
                        for (Vertex* v = first[bxj * Ny + byj]; v; v = v->next) // in first we get the starting vertex from each box
                            if (v != &polygons[p].vertices[k] && v != &polygons[p].vertices[i] && v != &polygons[p].vertices[j])
                                interaction(v, &polygons[p].vertices[i], &polygons[p].vertices[k], &polygons[p].vertices[j]);
            }
        }
        
        // time integration (semi-implicit Euler method)
        #pragma omp parallel for
        for (std::size_t p = Nr; p < polygons.size(); ++p) {
            // move vertices
            for (auto& v : polygons[p].vertices) {
                v.v.add(dt, v.a); // update vertex velocity
                v.r.add(dt, v.v); // update vertex position
            }
            // apply area growth rate
            polygons[p].area(); // compute the new polygon area
            double alpha_new = growthRate(polygons[p], polygons[p].c, polygons[p].grad_c, t);
            if (polygons[p].A > beta * polygons[p].A0 || alpha_new < 0) {
                polygons[p].A0 += alpha_new * dt;
            }
        }
        t += dt; // advance the time
    }
    
    // place all vertices into boxes
    void boxes() {
        // compute the global bounding box and maximum squared edge length
        double xmin = polygons[0].vertices[0].r.x, xmax = xmin;
        double ymin = polygons[0].vertices[0].r.y, ymax = ymin;
        double l2max = 0;
        #pragma omp parallel for reduction(min:xmin,ymin) reduction(max:xmax,ymax,l2max)
        for (std::size_t p = 0; p < polygons.size(); ++p) {
            for (std::size_t i = polygons[p].vertices.size() - 1, j = 0; j < polygons[p].vertices.size(); i = j++) {
                if      (polygons[p].vertices[i].r.x < xmin) xmin = polygons[p].vertices[i].r.x;
                else if (polygons[p].vertices[i].r.x > xmax) xmax = polygons[p].vertices[i].r.x;
                if      (polygons[p].vertices[i].r.y < ymin) ymin = polygons[p].vertices[i].r.y;
                else if (polygons[p].vertices[i].r.y > ymax) ymax = polygons[p].vertices[i].r.y;
                l2max = std::max(l2max, (polygons[p].vertices[j].r - polygons[p].vertices[i].r).length2());
            }
        }
        bbox.x0 = xmin;
        bbox.y0 = ymin;
        bbox.x1 = xmax;
        bbox.y1 = ymax;
        
        // place vertices in boxes
        bs = std::sqrt(l2max) + drmax; // box size
        Nx = (xmax - bbox.x0) / bs + 3; // number of boxes in x direction (with an extra column on both ends)
        Ny = (ymax - bbox.y0) / bs + 3; // number of boxes in y direction (with an extra row on both ends)
        first.assign(Nx * Ny, 0); // clear the boxes
        for (auto& p : polygons) {
            for (auto& v : p.vertices) {
                const std::size_t bx = (v.r.x - bbox.x0) / bs + 1; // box index in x direction
                const std::size_t by = (v.r.y - bbox.y0) / bs + 1; // box index in y direction
                const std::size_t b = bx * Ny + by; // global box index
                v.next = first[b];
                first[b] = &v; // place this vertex in front of the list in this box
            }
        }
    }
    
    // polygon-polygon interaction
    void interaction(Vertex* v0, Vertex* v1, Vertex* v1n0, Vertex* v1n1) {
        const Polygon& p0 = polygons[v0->p]; // polygon 0
        const Polygon& p1 = polygons[v1->p]; // polygon 1
        // select the closer of the two edges
        Vertex* v1n [2] = {v1n0, v1n1}; // pointers to neighbors of vertex 1
        double xi [2], xit [2], dr2 [2];
        Point dr [2];
        for (unsigned int j = 0; j < 2; ++j)
            dr2[j] = point_edge_dist2(v0->r, v1->r, v1n[j]->r, xi[j], xit[j], dr[j]);
        const unsigned int j = dr2[0] > dr2[1];
        
        if (xi[j] <= 1 && dr2[j] < drmax * drmax) {
            Vertex* v2 = v1n[j]; // pointer to vertex 2
            
            // check if arclength distance is large enough to allow for self-contact
            if (lmin < h && v0->p == v1->p) {
                auto& v = p1.vertices;
                const std::size_t j1 = v1 - v.data(), j2 = v2 - v.data();
                for (unsigned int d = 0; d < 2; ++d) {
                    double l = 0;
                    for (std::size_t i = v0 - v.data(), k; l <= h && i != j1 && i != j2; i = k) {
                        k = (i + v.size() + 2*d-1) % v.size();
                        l += (v[i].r - v[k].r).length();
                    }
                    if (l <= h) return;
                }
            }
            
            // trilinear traction-separation law including viscous damping
            const double dr_abs = std::sqrt(dr2[j]);
            const Point dv = v0->v - v1->v - xit[j] * (v2->v - v1->v);
            const double dndv_dr = dr[j] * dv / dr2[j];
            double da_dr = -cc * dndv_dr; // viscous dashpot
            if (dr_abs < h) {
                da_dr += kr * (h / dr_abs - 1); // linear repulsion
            }
            else if (v0->p != v1->p) {
                const double kh_new = adhesionStiffness(p0, p1, p0.c, p1.c, p0.grad_c, p1.grad_c, t);
                if (dr_abs < h + sh)
                    da_dr += kh_new * (h / dr_abs - 1); // linear adhesion
                else
                    da_dr += kh_new * sh / ss * (1 - drmax / dr_abs); // softening adhesion
            }
            Point a = da_dr * dr[j]; // acceleration vector
            
            // add dynamic Coulomb friction
            const Point dvt = dv - dndv_dr * dr[j];
            const double dvt_abs = dvt.length();
            const double mu_new = dynamicFriction(p0, p1, p0.c, p1.c, p0.grad_c, p1.grad_c, t);
            if (dvt_abs > 0 && v0->p != v1->p)
                a.add(-mu_new * std::abs(da_dr) * dr_abs / dvt_abs, dvt);
            
            // distribute the acceleration to the involved vertices
            v0->a.add(1, a);
            v1->a.add(xit[j]-1, a);
            v2->a.add(-xit[j], a);
        }
    }
    
    void output(const std::size_t f) {
        // print the frame number, simulated time, number of polygons, and total number of vertices
        std::size_t Nv = 0;
        for (auto& p : polygons)
            Nv += p.vertices.size();
        std::cout << "frame " << f << ", t=" << t << ", " << polygons.size() << " polygons, " << Nv << " vertices\n";
        
        // write a VTK polygon file containing the areas, perimeters, and coordination numbers
        char name [16];
        snprintf(name, 16, "frame%06zu.vtp", f);
        std::ofstream file(output_folder + "/" + name);
        if (!file) {
            std::cerr << "Error: Could not open file " << output_folder << "/" << name << std::endl;
            return;
        }
        file << "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
        file << "    <PolyData>\n";
        file << "        <Piece NumberOfPoints=\"" << Nv << "\" NumberOfPolys=\"" << polygons.size() << "\">\n";
        file << "            <CellData>\n";
        file << "                <DataArray type=\"Float64\" Name=\"area\" format=\"ascii\">\n";
        for (auto& p : polygons)
            file << p.A << " ";
        file << "\n";
        file << "                </DataArray>\n";
        if (Output::c) {
            for (std::size_t i = 0; i < NUM_SPECIES; ++i) {
                file << "                <DataArray type=\"Float64\" Name=\"c" << i+1 << "\" format=\"ascii\">\n";
                for (auto& p : polygons)
                    file << p.c[i] << " ";
                file << "\n";
                file << "                </DataArray>\n";
            }
        }
        if (Output::grad_c) {
            for (std::size_t i = 0; i < NUM_SPECIES; ++i) {
                file << "                <DataArray type=\"Float64\" Name=\"grad_c" << i+1 << "\" NumberOfComponents=\"2\" format=\"ascii\">\n";
                for (auto& p : polygons)
                    file << p.grad_c[i].x << " " << p.grad_c[i].y << " ";
                file << "\n";
                file << "                </DataArray>\n";
            }
        }
        if (Output::D) {
            for (std::size_t i = 0; i < NUM_SPECIES; ++i) {
                file << "                <DataArray type=\"Float64\" Name=\"D" << i+1 << "\" format=\"ascii\">\n";
                for (auto& p : polygons)
                    file << p.D[i] << " ";
                file << "\n";
                file << "                </DataArray>\n";
            }
        }
        if (Output::k) {
            for (std::size_t i = 0; i < NUM_KINETIC; ++i) {
                file << "                <DataArray type=\"Float64\" Name=\"k" << i+1 << "\" format=\"ascii\">\n";
                for (auto& p : polygons)
                      file << p.k[i] << " ";
                file << "\n";
                file << "                </DataArray>\n";
            }
        }
        if (Output::cell_type) {
            file << "                <DataArray type=\"Int64\" Name=\"cell_type\" format=\"ascii\">\n";
            for (auto& p : polygons)
                file << p.cell_type << " ";
            file << "\n";
            file << "                </DataArray>\n";
        }
        file << "                <DataArray type=\"Float64\" Name=\"perimeter\" format=\"ascii\">\n";
        for (auto& p : polygons)
            file << p.perimeter() << " ";
        file << "\n";
        file << "                </DataArray>\n";
        file << "                <DataArray type=\"Float64\" Name=\"neighbors\" format=\"ascii\">\n";
        boxes(); // place all vertices into boxes
        #pragma omp parallel for ordered schedule(static,1)
        for (std::size_t p = 0; p < polygons.size(); ++p) {
            double xi, xit;
            Point dr;
            std::set<std::size_t> n; // set of neighboring polygons
            for (std::size_t i = polygons[p].vertices.size() - 1, k = i - 1, j = 0; j < polygons[p].vertices.size(); k = i, i = j++) {
                const std::size_t bxi = (polygons[p].vertices[i].r.x - bbox.x0) / bs + 1; // box index in x direction
                const std::size_t byi = (polygons[p].vertices[i].r.y - bbox.y0) / bs + 1; // box index in y direction
                for (std::size_t bxj = bxi - 1; bxj <= bxi + 1; ++bxj)
                    for (std::size_t byj = byi - 1; byj <= byi + 1; ++byj)
                        for (Vertex* v = first[bxj * Ny + byj]; v; v = v->next)
                            if (p != v->p) // do not count contact of a polygon with itself
                                if (point_edge_dist2(v->r, polygons[p].vertices[i].r, polygons[p].vertices[k].r, xi, xit, dr) < drmax * drmax ||
                                    point_edge_dist2(v->r, polygons[p].vertices[i].r, polygons[p].vertices[j].r, xi, xit, dr) < drmax * drmax)
                                    n.insert(v->p);
            }
            #pragma omp ordered
            file << n.size() << " ";
        }
        file << "\n";
        file << "                </DataArray>\n";
        file << "            </CellData>\n";
        file << "            <Points>\n";
        file << "                <DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        for (auto& p : polygons)
            for (auto& v : p.vertices)
                file << v.r.x << " " << v.r.y << " " << p.phase << " ";
        file << "\n";
        file << "                </DataArray>\n";
        file << "            </Points>\n";
        file << "            <Polys>\n";
        file << "                <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n";
        for (std::size_t i = 0; i < Nv; ++i)
            file << i << " ";
        file << "\n";
        file << "                </DataArray>\n";
        file << "                <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n";
        std::size_t offset = 0;
        for (auto& p : polygons)
            file << (offset += p.vertices.size()) << " ";
        file << "\n";
        file << "                </DataArray>\n";
        file << "            </Polys>\n";
        file << "        </Piece>\n";
        file << "    </PolyData>\n";
        file << "</VTKFile>\n";
    }

    // write the current ensemble state to file for later entry point
    void write_OFF(std::string filename) const {
        std::ofstream file(output_folder + "/" + filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << "/" << filename << std::endl;
            return;
        }
        
        // count all vertices
        std::size_t Nv = 0;
        for (auto& p : polygons)
            Nv += p.vertices.size();
        
        file << "OFF\n";
        file << Nv << " " << polygons.size() << " 0\n";
        
        // write all vertices (ensure no duplicate vertices)
        std::vector<Point> allPoints;
        for (const auto& p : polygons)
            for (const auto& v : p.vertices)
                allPoints.push_back(v.r);
        
        // Removing duplicate points and mapping original indices to new indices
        std::vector<Point> uniquePoints;
        std::map<Point, int> pointIndexMap;
        std::size_t index = 0;
        for (const auto& point : allPoints) {
            if (pointIndexMap.find(point) == pointIndexMap.end()) {
                uniquePoints.push_back(point);
                pointIndexMap[point] = index++;
            }
        }
        
        // write unique points to file
        for (const auto& point : uniquePoints)
            file << point.x << " " << point.y << " 0\n";
        
        // write polygons using the indices of unique points
        for (const auto& p : polygons) {
            file << p.vertices.size();
            for (const auto& v : p.vertices)
                file << " " << pointIndexMap[v.r];
            file << "\n";
        }
        file.close();
    }
};

#endif
