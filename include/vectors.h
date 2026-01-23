#ifndef VECTORS_H
#define VECTORS_H

#include <stddef.h>

double normalize_vector(double *v, size_t n);

double sum_vector(double *v, size_t n);

double dot_product(double *a, double *b, size_t n);

double standard_deviation(size_t N, const double *x);

double max_vector(double *v, size_t n);

double min_vector(double *v, size_t n);

#endif
