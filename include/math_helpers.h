// math_utils.h
#ifndef MATH_HELPERS_H
#define MATH_HELPERS_H

#include <math.h>
#include <stddef.h>
#include <string.h>

size_t *ilogspace_unique(double log10_ini, double log10_end, size_t n,
                         size_t *out_n_unique);

static inline int sgn(double x) {

    if (x > 0)
        return 1;
    if (x < 0)
        return -1;
    return 0;
}

double area_trapezoidal(const double *x, const double *y, size_t n);

static inline void welford_update(const double *h, double *mean, double *M2,
                                  size_t dof, size_t *count) {
    (*count)++;
    for (size_t i = 0; i < dof; i++) {

        // --- Update for E ---
        double e = h[i];
        double d1 = e - mean[i];
        mean[i] += d1 / (double)*count;
        double d2 = e - mean[i];
        M2[i] += d1 * d2;
    }
}

static inline void welford_reset(double *mean, double *M2, size_t n) {
    memset(mean, 0, n * sizeof *mean);
    memset(M2, 0, n * sizeof *M2);
}

static inline double welford_std(const double *M2, size_t dof, size_t count) {

    double sigma = 0.0;

    if (count == 0)
        return -1;

    for (size_t i = 0; i < dof; i++) {
        sigma += sqrt(M2[i] / count);
    }
    sigma /= dof;

    return sigma;
}

#endif
