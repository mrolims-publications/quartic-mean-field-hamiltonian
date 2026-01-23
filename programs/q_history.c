#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/math_helpers.h"
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

    if (argc != 5) {
        fprintf(
            stderr,
            "Usage: %s <dof> <eps> <total_time> <seed>\n\n"
            "Arguments:\n"
            "  <dof>         Number of degrees of freedom (positive integer)\n"
            "  <eps>         Specific energy ε = E / N\n"
            "  <total_time>  Total integration time\n"
            "  <seed>        Random number generator seed\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    size_t dof = get_size_t(argc, argv, 1);
    double eps = get_double(argc, argv, 2);
    double total_time = get_double(argc, argv, 3);
    unsigned int seed = get_uint(argc, argv, 4);

    double p_range = 7.5 * sqrt(2 * eps);

    size_t n_steps = (size_t)(total_time / TIME_STEP);
    size_t sample_ini = (size_t)(1e4 / TIME_STEP);

    double *q = xmalloc(dof, sizeof *q);
    double *p = xmalloc(dof, sizeof *p);
    generate_waterbag_initial_condition(dof, eps, q, p, seed);

    char filename[512];
    snprintf(filename, sizeof(filename),
             "data/q_history_dof=%zu_eps=%.6f_T=%.1f_seed=%d.dat", dof, eps,
             total_time, seed);
    FILE *fp = fopen(filename, "w");

    size_t n_unique;
    size_t *sample_steps =
        ilogspace_unique(log10(sample_ini), log10(n_steps), 5000, &n_unique);

    double *dV = xmalloc(dof, sizeof *dV);
    double *pdf = xcalloc(N_BINS, sizeof *pdf);
    double *pdf_norm = xcalloc(N_BINS, sizeof *pdf_norm);
    double *X = xmalloc(1000, sizeof *X);
    double *Y = xmalloc(1000, sizeof *Y);
    double bin_width = 2 * p_range / N_BINS;
    size_t sample_idx = 0;
    double A0 = 0.9, beta0 = 1.0, q0 = 1.1;
    for (size_t n = 0; n < n_steps; n++) {
        yoshida4_step(q, p, dof, dV);

        if ((n + 1) % SAMPLE_STRIDE == 0) {
            for (size_t j = 0; j < dof; j++) {
                double pj = p[j];
                if (pj < -p_range || pj > p_range) {
                    fprintf(stderr,
                            "ERROR: p = %.5f is outside p_range for dof = %zu "
                            "and eps = %.1f\n",
                            pj, dof, eps);
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

        if (n + 1 == sample_steps[sample_idx]) {
            double total_counts = 0.0;
            size_t positive_counts = 0;
            for (size_t i = 0; i < N_BINS; i++) {
                total_counts += pdf[i];
                if (pdf[i] > 0)
                    positive_counts++;
            }

            for (size_t i = 0; i < N_BINS; i++)
                pdf_norm[i] = pdf[i] / (total_counts * bin_width);

            double *x_fit = xmalloc(positive_counts, sizeof *x_fit);
            double *y_fit = xmalloc(positive_counts, sizeof *y_fit);
            double *weights = xmalloc(positive_counts, sizeof *weights);

            size_t zero_bin = (size_t)((0.0 + p_range) / bin_width);
            double P0 = (zero_bin < N_BINS) ? pdf_norm[zero_bin] : 1.0;
            size_t count = 0;
            for (size_t i = 0; i < N_BINS; i++) {
                double bin_center = -p_range + (i + 0.5) * bin_width;
                if (pdf_norm[i] > 0) {
                    x_fit[count] = bin_center * P0;
                    y_fit[count] = (pdf_norm[i] / P0);
                    weights[count] = pdf[i];
                    count++;
                }
            }

            QGaussianFitResult res =
                fit_q_gaussian(x_fit, y_fit, count, A0, beta0, q0);

            FitMetrics metrics = compute_fit_metrics(x_fit, y_fit, res, count);

            double area_num = area_trapezoidal(x_fit, y_fit, count);

            double xmin = x_fit[0];
            double xmax = x_fit[count - 1];

            for (size_t i = 0; i < 1000; i++) {
                X[i] = xmin + i * (xmax - xmin) / (1000 - 1);
                Y[i] = q_gaussian(X[i], res.A, res.beta, res.q);
            }
            double area_fit = area_trapezoidal(X, Y, 1000);

            fprintf(fp, "%.16f %.16f %.16e %.16e %.16e %.16e %.16e\n",
                    (n + 1) * TIME_STEP, res.q, res.q_err, metrics.mse,
                    metrics.mae, area_num, area_fit);

            free(x_fit);
            free(y_fit);
            free(weights);
            sample_idx++;
            A0 = res.A;
            beta0 = res.beta;
            q0 = res.q;
        }
    }

    free(dV);
    free(q);
    free(p);
    free(sample_steps);
    free(pdf);
    free(pdf_norm);
    free(X);
    free(Y);

    fclose(fp);
    return 0;
}
