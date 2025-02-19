// PolyMorph
// Copyright (c) 2024-2025
// Nicolas Müller, nicolmueller@ethz.ch
// Roman Vetter, vetterro@ethz.ch
// ETH Zurich

#ifndef DOMAIN_H
#define DOMAIN_H

#include "geometry.h"
#include "param.h"

struct Domain {
    double x0, y0, x1, y1;

    Domain() : x0(0), y0(0), x1(1), y1(1) {}
    Domain(double x0, double y0, double x1, double y1) : 
            x0(x0), y0(y0), x1(x1), y1(y1) {}

    double width() const { return x1 - x0; }
    double height() const { return y1 - y0; }

    bool contains(const Point& r) const {
        return r.x >= x0 && r.x <= x1 && r.y >= y0 && r.y <= y1;
    }
};

#endif // DOMAIN_H