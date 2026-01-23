#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/model.h"
#include "../include/parse_args.h"
#include "../include/q_gaussian_fit.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_BINS 1000
#define SAMPLE_STRIDE 5

int main(int argc, char *argv[]) {

    size_t dof = get_size_t(argc, argv, 1);
    double eps = get_double(argc, argv, 2);
    double total_time = get_double(argc, argv, 3);
    unsigned int seed = get_uint(argc, argv, 4);

    double p_range = 7.5 * sqrt(2 * eps);

    size_t n_steps = (size_t)(total_time / TIME_STEP);

    double *q = xmalloc(dof, sizeof *q);
    double *p = xmalloc(dof, sizeof *p);
    generate_waterbag_initial_condition(dof, eps, q, p, seed);

    char filename[512];
    snprintf(filename, sizeof(filename),
             "data/p_PDF_dof=%zu_eps=%.6f_T=%.1f_seed=%d.dat", dof, eps,
             total_time, seed);
    FILE *fp = fopen(filename, "w");

    double *dV = xmalloc(dof, sizeof *dV);
    double *pdf = xcalloc(N_BINS, sizeof *pdf);
    double bin_width = 2 * p_range / N_BINS;
    for (size_t n = 0; n < n_steps; n++) {
        yoshida4_step(q, p, dof, dV);

        if ((n + 1) % SAMPLE_STRIDE == 0) {
            for (size_t j = 0; j < dof; j++) {
                double pj = p[j];
                if (pj < -p_range || pj > p_range) {
                    fprintf(stderr, "ERROR: p = %.5f is outside p_range\n", pj);
                    exit(EXIT_FAILURE);
                }
                int bin = (int)((pj + p_range) / bin_width);
                if (bin < 0)
                    bin = 0;
                if ((size_t)bin >= N_BINS)
                    bin = N_BINS - 1;
                pdf[bin]++;
            }
        }
    }

    double total_counts = 0.0;
    size_t positive_counts = 0.0;
    for (size_t i = 0; i < N_BINS; i++) {
        total_counts += pdf[i];
        if (pdf[i] > 0)
            positive_counts++;
    }

    for (size_t i = 0; i < N_BINS; i++)
        pdf[i] /= (total_counts * bin_width);

    size_t zero_bin = (size_t)((0.0 + p_range) / bin_width);
    double P0 = (zero_bin < N_BINS) ? pdf[zero_bin] : 1.0;

    double *x_fit = xmalloc(positive_counts, sizeof *x_fit);
    double *y_fit = xmalloc(positive_counts, sizeof *y_fit);
    size_t count = 0;
    for (size_t i = 0; i < N_BINS; i++) {
        double bin_center = -p_range + (i + 0.5) * bin_width;
        double x_value = bin_center * P0;                // p * P(0)
        double y_value = (P0 > 0) ? (pdf[i] / P0) : 0.0; // P(p) / P(0)
        fprintf(fp, "%.16f %.16f\n", x_value, y_value);

        if (y_value > 0) {
            x_fit[count] = x_value;
            y_fit[count] = y_value;
            count++;
        }
    }

    const double A0 = 0.9, beta0 = 1.0, q0 = 1.1;
    QGaussianFitResult res = fit_q_gaussian(x_fit, y_fit, count, A0, beta0, q0);

    FitMetrics metrics = compute_fit_metrics(x_fit, y_fit, res, count);

    printf("q = %.16f\n", res.q);
    printf("std_q = %.16f\n", res.q_err);
    printf("mse = %.16f\n", metrics.mse);
    printf("mae = %.16f\n", metrics.mae);

    free(dV);
    free(q);
    free(p);
    free(pdf);

    return 0;
}
