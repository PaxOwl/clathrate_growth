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