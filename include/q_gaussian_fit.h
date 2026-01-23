#ifndef Q_GAUSSIAN_FIT_H
#define Q_GAUSSIAN_FIT_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    double A;     /* fitted amplitude */
    double beta;  /* fitted beta */
    double q;     /* fitted q */
    double q_err; /* std error of q (sqrt of covariance diag) */
    int status;   /* GSL status code (0 == success) */
} QGaussianFitResult;

typedef struct {
    double mse;
    double mae;
} FitMetrics;

double q_gaussian(double p, double A, double beta, double q);

QGaussianFitResult fit_q_gaussian(const double *p, const double *P, size_t n,
                                  double A0, double beta0, double q0);

FitMetrics compute_fit_metrics(double *x, double *y, QGaussianFitResult res,
                               size_t n);

#ifdef __cplusplus
}
#endif

#endif /* Q_GAUSSIAN_FIT_H */
