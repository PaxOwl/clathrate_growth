#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"


void periodic_conditions(double* delta, double* box, void* out) {
    double* outdata = (double*) out;
    for(int i; i < 3; i++) {
        delta[i] = delta[i] - round(delta[i] / box[i]) * box[i];
    }
}

double* distance(double* p1, double* p2, double* box, bool periodic, double* results) {
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    double delta[3] = {dx, dy, dz};

    if (periodic == true) {
        results = periodic_conditions(delta, box);
    }
}
