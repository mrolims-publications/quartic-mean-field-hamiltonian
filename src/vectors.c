#include "../include/vectors.h"
#include <float.h>
#include <math.h>

double normalize_vector(double *v, size_t n) {
    double norm = 0.0;
    for (size_t i = 0; i < n; i++)
        norm += v[i] * v[i];
    norm = sqrt(norm);
    for (size_t i = 0; i < n; i++)
        v[i] /= norm;
    return norm;
}

double sum_vector(double *v, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
        sum += v[i];
    return sum;
}

double dot_product(double *a, double *b, size_t n) {
    double sum = 0.0;
    for (size_t i = 0; i < n; i++)
        sum += a[i] * b[i];

    return sqrt(sum);
}

double standard_deviation(size_t N, const double *x) {

    double mean = 0.0;
    double std = 0.0;

    // Compute means
    for (size_t i = 0; i < N; i++) {
        mean += x[i];
    }
    mean /= N;

    // Compute variance
    for (size_t i = 0; i < N; i++) {
        double diff = x[i] - mean;
        std += diff * diff;
    }
    std = sqrt(std / N);

    return std;
}

double max_vector(double *v, size_t n) {

    double max = -DBL_MAX;
    for (size_t i = 0; i < n; i++)
        max = fmax(v[i], max);

    return max;
}

double min_vector(double *v, size_t n) {

    double min = DBL_MAX;
    for (size_t i = 0; i < n; i++)
        min = fmin(v[i], min);

    return min;
}
