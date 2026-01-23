#include "../include/math_helpers.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

size_t *ilogspace_unique(double log10_ini, double log10_end, size_t n,
                         size_t *out_n_unique) {
    size_t *tmp = malloc(n * sizeof *tmp);
    if (!tmp) {
        fprintf(stderr, "malloc has failed\n");
        exit(EXIT_FAILURE);
    }

    size_t n_unique = 0;
    size_t last_val = (size_t)-1;

    for (size_t i = 0; i < n; i++) {
        double exponent = log10_ini + i * (log10_end - log10_ini) / (n - 1);
        size_t val = (size_t)pow(10.0, exponent);
        if (val != last_val) {
            tmp[n_unique++] = val;
            last_val = val;
        }
    }

    size_t *result = malloc(n_unique * sizeof *result);
    if (!result) {
        fprintf(stderr, "malloc has failed\n");
        exit(EXIT_FAILURE);
    }

    for (size_t i = 0; i < n_unique; i++)
        result[i] = tmp[i];
    free(tmp);

    *out_n_unique = n_unique;
    return result;
}

double area_trapezoidal(const double *x, const double *y, size_t n) {
    if (n < 2)
        return 0.0; // not enough points

    double area = 0.0;

    for (size_t i = 0; i < n - 1; i++) {
        double dx = x[i + 1] - x[i];
        double avg = 0.5 * (y[i] + y[i + 1]);
        area += avg * dx;
    }

    return area;
}
