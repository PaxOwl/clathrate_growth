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
    double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    theta[0] = acos(dot);
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
                        double l_lim, double h_lim,
                        int n_size, double *vec, long *out) {
    double temp_dst;
    size_t counter = 0;
    size_t i;
    for (i = 0; i < n_size; i++) {
        if (center[0] == neighbours[i][0]
            && center[1] == neighbours[i][1]
            && center[2] == neighbours[i][2]) {
            continue;
        }
        distance(center, neighbours[i], box, vec);
        temp_dst = sqrt(pow(vec[0], 2) + pow(vec[1], 2) + pow(vec[2], 2));
        if (temp_dst <= h_lim && temp_dst >= l_lim) {
            out[counter] = i;
            counter++;
        }
    }
}

void neighbours(double (*centers)[3], double (*neighbours)[3], double *box,
                double l_lim, double h_lim, int c_size, int n_size,
                double *vec, long (*out)[n_size]) {
    size_t i;
    for (i = 0; i < c_size; i++) {
        printf("%lu\n", i);
        nearest_neighbours(centers[i], neighbours, box,
                           h_lim, l_lim, n_size, vec, out[i]);
    }
}

void aop(double *center, double (*neighbours)[3], double *box,
         int n_size, double *vec1, double *vec2, double *theta, double *angles){
    size_t i;
    size_t n_max = n_size - 1;
    size_t n_min = 0;
    size_t iter = 0;
    do {
        for (i = n_min; i < n_max; i++) {
            distance(center, neighbours[n_min], box, vec1);
            distance(center, neighbours[i + 1], box, vec2);
            angle(vec1, vec2, theta);
            angles[i - n_min + n_min * n_size] = theta[0];
        }
        n_min++;
    } while (n_min < n_max);
}