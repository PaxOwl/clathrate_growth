#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"


void periodic_conditions(double* delta, double* box) {
    for(int i; i < 3; i++) {
        delta[i] = delta[i] - round(delta[i] / box[i]) * box[i];
    }
}

void distance(double* p1, double* p2, double* box, double* out) {
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    out[0] = dx;
    out[1] = dy;
    out[2] = dz;

    periodic_conditions(out, box);
}

void norm_vec(double* vec) {
    double norm = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
    for (int i; i < 3; i++) {
        vec[i] /= norm;
    }
}

void angle(double* v1, double* v2, double* theta) {
    double pi = 3.14159265359;
    double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    theta[0] = acos(dot) * 180 / pi;
}

void closest_atom(double* names, double* oxygens, double* hydrogens,
                  double* box, double* out_name) {
    double do1h1[3] = {0., 0., 0.};
    double do2h1[3] = {0., 0., 0.};
    double do1h2[3] = {0., 0., 0.};
    double do2h2[3] = {0., 0., 0.};
    double dsth1;
    double dsth2;

    distance(&oxygens[0], &hydrogens[0], box, do1h1);
    distance(&oxygens[1], &hydrogens[0], box, do2h1);
    distance(&oxygens[0], &hydrogens[1], box, do1h2);
    distance(&oxygens[1], &hydrogens[1], box, do2h2);

    dsth1 = sqrt(pow(do1h1[0], 2) + pow(do1h1[1], 2) + pow(do1h1[2], 2))
            + sqrt(pow(do1h1[0], 2) + pow(do2h1[1], 2) + pow(do2h1[2], 2));
    dsth2 = sqrt(pow(do1h2[0], 2) + pow(do1h2[1], 2) + pow(do1h2[2], 2))
            + sqrt(pow(do2h2[0], 2) + pow(do2h2[1], 2) + pow(do2h2[2], 2));

    if (dsth2 > dsth1) {
        out_name[0] = names[0];
    }
    else {
        out_name[0] = names[1];
    }
}
