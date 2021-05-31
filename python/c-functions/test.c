#include <stdio.h>


void test(int size, double (*matrice)[size]) {
    printf("%d\n", size);
    size_t i, j;
    for (i = 0; i < size; i++) {
        printf("loop 1 - %lu\n", i);
        for (j = 0; j < size; j++) {
            printf("loop 2 - %lu\n", j);
            printf("%f\n", matrice[i][j]);
        }
    }
}