#include <stdbool.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"


void periodic_conditions(double *delta, double *box) {
    size_t i;
    for(i = 0; i < 3; i++) {
        delta[i] = delta[i] - round(delta[i] / box[i]) * box[i];
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

void angle(double *v1, double *v2, double *theta) {
    double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    double norm1 = sqrt(pow(v1[0], 2) + pow(v1[1], 2) + pow(v1[2], 2));
    double norm2 = sqrt(pow(v2[0], 2) + pow(v2[1], 2) + pow(v2[2], 2));
    theta[0] = acos(dot / (norm1 * norm2));
}

void closest_atom(double *center, double *oxygens, double (*hydrogens)[3],
                  double *box, double *closest) {
    double do1h1[3] = {0., 0., 0.};
    double do2h1[3] = {0., 0., 0.};
    double do1h2[3] = {0., 0., 0.};
    double do2h2[3] = {0., 0., 0.};
    double dsth1;
    double dsth2;

    distance(center, hydrogens[0], box, do1h1);
    distance(oxygens, hydrogens[0], box, do2h1);
    distance(center, hydrogens[1], box, do1h2);
    distance(oxygens, hydrogens[1], box, do2h2);

    dsth1 = sqrt(pow(do1h1[0], 2) + pow(do1h1[1], 2) + pow(do1h1[2], 2))
            + sqrt(pow(do2h1[0], 2) + pow(do2h1[1], 2) + pow(do2h1[2], 2));
    dsth2 = sqrt(pow(do1h2[0], 2) + pow(do1h2[1], 2) + pow(do1h2[2], 2))
            + sqrt(pow(do2h2[0], 2) + pow(do2h2[1], 2) + pow(do2h2[2], 2));

    if (dsth2 >= dsth1) {
        closest[0] = hydrogens[0][0];
        closest[1] = hydrogens[0][1];
        closest[2] = hydrogens[0][2];
    }
    else {
        closest[0] = hydrogens[1][0];
        closest[1] = hydrogens[1][1];
        closest[2] = hydrogens[1][2];
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
        nearest_neighbours(centers[i], neighbours, box,
                           h_lim, l_lim, n_size, vec, out[i]);
    }
}

void aop(double *center, double (*neighbours)[3], double *box,
         int n_size, double *vec1, double *vec2, double *theta,
         double *angles, double *aop){
    double pi = 3.14159265359;
    size_t i;
    size_t iter = 0;
    size_t n_min = 0;
    do {
        for (i = n_min; i < n_size - 1; i++) {
            distance(center, neighbours[n_min], box, vec1);
            distance(center, neighbours[i + 1], box, vec2);
            angle(vec1, vec2, theta);
            angles[iter] = theta[0];
            iter++;
        }
        n_min++;
    } while (n_min < n_size + 1);
    for (i = 0; i < iter; i++) {
        aop[0] = aop[0] + pow(fabs(cos(angles[i])) * cos(angles[i])
                              + pow(cos(109.47 * pi / 180), 2), 2);
    }
}

void hydrogen_bonds(double *center, double (*oxygens)[3],
                    double (*hydrogens)[3], double *box,
                    double *vec1, double *vec2, double* vec3,
                    int ox_size, double *theta,
                    double *closest, double *bonds) {
    double pi = 3.14159265359;
    size_t i;
    double dst;
    for (i = 0; i < ox_size; i++) {
        if (center[0] == oxygens[i][0]
            && center[1] == oxygens[i][1]
            && center[2] == oxygens[i][2]) {
            continue;
        }
        distance(center, oxygens[i], box, vec1);
        dst = sqrt(pow(vec1[0], 2) + pow(vec1[1], 2) + pow(vec1[2], 2));
        if (dst > 0.25 && dst < 0.35) {
            closest_atom(center, oxygens[i], hydrogens, box, closest);
            distance(closest, center, box, vec2);
            distance(closest, oxygens[i], box, vec3);
            angle(vec2, vec3, theta);
            theta[0] = theta[0] * 180 / pi;
            if (theta[0] > 90 && theta[0] < 180) {
                bonds[i] = 1;
            }
        }
            theta[0] = 0;
    }
}

void clath_size(double *small, double *large, int n_small, int n_large,
          double *xmin, double *xmax) {
    for (size_t i; i < n_small; i++) {
        if (small[i] > xmax[0]) {
            xmax[0] = small[i];
        }
        if (small[i] < xmin[0]) {
            xmin[0] = small[i];
        }
    }
    for (size_t i; i < n_large; i++) {
        if (small[i] > xmax[0]) {
            xmax[0] = large[i];
        }
        if (small[i] < xmin[0]) {
            xmin[0] = large[i];
        }
    }

}
