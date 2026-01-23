#include "../include/q_gaussian_fit.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_multifit_nlinear.h>
#include <gsl/gsl_version.h>
#include <math.h>
#include <string.h>

/* Model and helpers */
static inline double inside_val(double p, double beta, double q) {
    return 1.0 - (1.0 - q) * beta * p * p;
}

double q_gaussian(double p, double A, double beta, double q) {
    const double t = inside_val(p, beta, q);
    if (t <= 0.0)
        return 0.0;
    return A * pow(t, 1.0 / (1.0 - q));
}

/* Data wrapper */
typedef struct {
    size_t n;
    const double *p;
    const double *P;
} data_t;

/* Residuals r_i = f(p_i; x) - P_i */
static int f_resid(const gsl_vector *x, void *params, gsl_vector *f) {
    const double A = gsl_vector_get(x, 0);
    const double beta = gsl_vector_get(x, 1);
    const double q = gsl_vector_get(x, 2);

    const data_t *d = (const data_t *)params;
    for (size_t i = 0; i < d->n; ++i) {
        const double yi = q_gaussian(d->p[i], A, beta, q);
        gsl_vector_set(f, i, yi - d->P[i]);
    }
    return GSL_SUCCESS;
}

/* Analytic Jacobian J_ij = d r_i / d x_j  (x = [A, beta, q]) */
static int f_jac(const gsl_vector *x, void *params, gsl_matrix *J) {
    const double A = gsl_vector_get(x, 0);
    const double beta = gsl_vector_get(x, 1);
    const double q = gsl_vector_get(x, 2);

    const data_t *d = (const data_t *)params;
    const double inv1mq = 1.0 / (1.0 - q);

    for (size_t i = 0; i < d->n; ++i) {
        const double pi = d->p[i];
        const double t = inside_val(pi, beta, q);

        double dA = 0.0, dbeta = 0.0, dq = 0.0;

        if (t > 0.0) {
            const double f = A * pow(t, inv1mq);
            const double t_pow_q_over_1mq =
                pow(t, (q)*inv1mq); /* t^{q/(1-q)} */
            /* ∂f/∂A = t^{1/(1-q)} */
            const double dfdA = pow(t, inv1mq);

            /* ∂f/∂β = -A * p^2 * t^{q/(1-q)} */
            const double dfdbeta = -A * (pi * pi) * t_pow_q_over_1mq;

            /* ∂f/∂q = f * [ (1/(1-q)^2) * ln t + (1/(1-q)) * (β p^2)/t ] */
            const double ln_t = log(t);
            const double term1 = (1.0 / ((1.0 - q) * (1.0 - q))) * ln_t;
            const double term2 = inv1mq * (beta * pi * pi) / t;
            const double dfdq = f * (term1 + term2);

            dA = dfdA;
            dbeta = dfdbeta;
            dq = dfdq;
        }
        /* residual r_i = f_i - P_i => Jacobian same as df (since dP_i/dx = 0)
         */
        gsl_matrix_set(J, i, 0, dA);
        gsl_matrix_set(J, i, 1, dbeta);
        gsl_matrix_set(J, i, 2, dq);
    }
    return GSL_SUCCESS;
}

QGaussianFitResult fit_q_gaussian(const double *p, const double *P, size_t n,
                                  double A0, double beta0, double q0) {
    QGaussianFitResult result;
    memset(&result, 0, sizeof result);

    if (!p || !P || n < 3) {
        result.status = GSL_EINVAL;
        return result;
    }

    const size_t n_params = 3;
    data_t d = (data_t){.n = n, .p = p, .P = P};

    gsl_multifit_nlinear_parameters params =
        gsl_multifit_nlinear_default_parameters();
    const gsl_multifit_nlinear_type *T = gsl_multifit_nlinear_trust;
    gsl_multifit_nlinear_workspace *w =
        gsl_multifit_nlinear_alloc(T, &params, n, n_params);
    if (!w) {
        result.status = GSL_ENOMEM;
        return result;
    }

    gsl_multifit_nlinear_fdf fdf;
    fdf.f = f_resid;
    fdf.df = f_jac; /* analytic Jacobian */
    fdf.fvv = NULL;
    fdf.n = n;
    fdf.p = n_params;
    fdf.params = &d;

    gsl_vector *x = gsl_vector_alloc(n_params);
    if (!x) {
        gsl_multifit_nlinear_free(w);
        result.status = GSL_ENOMEM;
        return result;
    }
    gsl_vector_set(x, 0, A0);
    gsl_vector_set(x, 1, beta0);
    gsl_vector_set(x, 2, q0);

    int status = gsl_multifit_nlinear_init(x, &fdf, w);
    if (status) {
        gsl_vector_free(x);
        gsl_multifit_nlinear_free(w);
        result.status = status;
        return result;
    }

    const double xtol = 1e-10;
    const double gtol = 1e-10;
    const double ftol = 0.0;

    size_t iter = 0;
    int info = 0;

    do {
        status = gsl_multifit_nlinear_iterate(w);
        if (status)
            break;

        status = gsl_multifit_nlinear_test(xtol, gtol, ftol, &info, w);
        ++iter;
    } while (status == GSL_CONTINUE && iter < 200);

    /* Fitted parameters */
    const gsl_vector *x_final = gsl_multifit_nlinear_position(w);
    result.A = gsl_vector_get(x_final, 0);
    result.beta = gsl_vector_get(x_final, 1);
    result.q = gsl_vector_get(x_final, 2);
    result.status = status;

    /* Covariance: cov = s^2 * (J^T J)^{-1}, with s^2 = chi2/(n-p) */
    gsl_matrix *covar = gsl_matrix_alloc(n_params, n_params);
    if (covar) {
        const gsl_matrix *J = gsl_multifit_nlinear_jac(w);
        gsl_multifit_nlinear_covar(J, 0.0, covar);

        const gsl_vector *r = gsl_multifit_nlinear_residual(w);
        double chisq = 0.0;
        gsl_blas_ddot(r, r, &chisq);

        const double dof = (double)n - (double)n_params;
        const double s2 = (dof > 0) ? chisq / dof : 0.0;

        /* scale covariance by s^2 */
        for (size_t i = 0; i < n_params; ++i)
            for (size_t j = 0; j < n_params; ++j)
                gsl_matrix_set(covar, i, j, gsl_matrix_get(covar, i, j) * s2);

        result.q_err = sqrt(fabs(gsl_matrix_get(covar, 2, 2)));
        gsl_matrix_free(covar);
    } else {
        result.q_err = NAN;
    }

    gsl_vector_free(x);
    gsl_multifit_nlinear_free(w);
    return result;
}

FitMetrics compute_fit_metrics(double *x, double *y, QGaussianFitResult res,
                               size_t n) {
    double mse = 0.0;
    double mae = 0.0;

    // Compute mean of y
    for (size_t i = 0; i < n; i++) {
        double y_fitted = q_gaussian(x[i], res.A, res.beta, res.q);
        double err = y[i] - y_fitted;

        mse += err * err;
        mae += fabs(err);
    }

    mse /= n;
    mae /= n;

    FitMetrics metrics;
    metrics.mse = mse;
    metrics.mae = mae;

    return metrics;
}
