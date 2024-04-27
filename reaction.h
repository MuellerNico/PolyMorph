#include <vector>
#include "const.h"


struct Reaction {
    std::vector<double> operator()(std::vector<double> u, std::vector<double> k) {
        std::vector<double> r(NUM_SPECIES);
        return r;
    }
};

struct LinearDegradation : Reaction {
    std::vector<double> operator()(std::vector<double> u, std::vector<double> k) {
        std::vector<double> r(NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; i++) {
            r[i] = -k[i]*u[i];
        }
        return r;
    }
};

struct Inhibition : Reaction {
    std::vector<double> operator()(std::vector<double> u, std::vector<double> k) {
        std::vector<double> r(NUM_SPECIES);
        r = {
            -k[0]*u[0] - k[2]*u[0]*u[1],
            -k[1]*u[1] - k[2]*u[0]*u[1]
        };
        return r;
    }
};