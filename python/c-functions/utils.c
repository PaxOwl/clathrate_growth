#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"


void periodic_conditions(double *delta, double *box) {
    size_t i;
    for(i = 0; i < 3; i++) {
        delta[i] = delta[i] - round(delta[i] / box[i])  *box[i];
    }
}

void distance(double *p1, double *p2, double *box, double *vec) {
    double dx = p2[0] - p1[0];
    double dy = p2[1] - p1[1];
    double dz = p2[2] - p1[2];
    printf("oui");
    vec[0] = dx;
    vec[1] = dy;
    vec[2] = dz;
    periodic_conditions(vec, box);
}

void norm_vec(double *vec) {
    double norm = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
    size_t i;
    for (i = 0; i < 3; i++) {
        vec[i] /= norm;
    }
}

void angle(double *v1, double *v2, double *theta) {
    double pi = 3.14159265359;
    double dot = v1[0]  *v2[0] + v1[1]  *v2[1] + v1[2]  *v2[2];
    theta[0] = acos(dot)  *180 / pi;
}

void closest_atom(double *names, double *oxygens, double *hydrogens,
                  double *box, double *out_name) {
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

void nearest_neighbours(double *center, double (*neighbours)[3], double *box,
                        double limit, int n_size, double *vec, double *out) {
    double dst;
    int counter = 0;
    size_t i;
    for (i = 0; i < n_size; i++) {
        distance(center, neighbours[i], box, vec);
        dst = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
        if (dst <= limit) {
            printf("Condition passed for i = %lu\n", i);
            out[counter] = i;
            counter++;
        }
    }
}

void neighbours(double (*centers)[3], double (*neighbours)[3], double *box,
                double limit, int c_size, int n_size,
                double *vec, double (*out)[c_size]) {
    size_t i;
    for (i = 0; i < c_size; i++) {
        nearest_neighbours(centers[i], neighbours, box,
                           limit, n_size, vec, out[i]);
    }
}