#include "../include/lyapunov.h"
#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/model.h"
#include "../include/vectors.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void LLE_history_stop_condition(double *q, double *p, size_t dof,
                                double total_time, size_t *sample_steps,
                                double *lle, double rtol, double atol,
                                size_t num_vals, FILE *fp) {

    size_t n_steps = (size_t)(total_time / TIME_STEP);
    size_t dim = 2 * dof;
    double *X = xmalloc(dim, sizeof *X);
    double *dV = xmalloc(dof, sizeof *dV);
    double *q2 = xmalloc(dof, sizeof *q2);
    double *Hdq = xmalloc(dof, sizeof *Hdq);

    for (size_t i = 0; i < dim; i++) {
        X[i] = (double)rand() / (double)RAND_MAX;
    }
    normalize_vector(X, dim);

    double *lle_vals = xmalloc(num_vals, sizeof *lle_vals);

    double H_ref = T(p, dof) + V(q, dof);

    double sum_log = 0.0;
    size_t sample_idx = 0;
    for (size_t n = 0; n < n_steps; n++) {
        double time = (n + 1) * TIME_STEP;
        yoshida4_step_traj_tan(q, p, q2, X, dV, Hdq, dof, 1);

        double norm = normalize_vector(X, dim);
        sum_log += log(norm);
        lle_vals[n % num_vals] = (sum_log / time);
        double mle_avg = 0.0;
        double mle_std = 0.0;
        if (n > num_vals) {
            mle_avg = sum_vector(lle_vals, num_vals) / num_vals;
            mle_std = standard_deviation(num_vals, lle_vals);
        }
        double standard_error = mle_std / sqrt((double)num_vals);

        if (standard_error > 0 &&
            (standard_error < atol ||
             (fabs(mle_avg) > 0.0 && standard_error / fabs(mle_avg) < rtol))) {
            printf("%.16f\n", standard_error);
            if (fp) {
                double E = T(p, dof) + V(q, dof);
                double Er = fabs(E - H_ref) / H_ref;
                fprintf(fp, "%.16f %.16e %.16e %.16e\n", time, Er,
                        lle_vals[n % num_vals], standard_error);
                fflush(fp);
            }
            break;
        }

        if (n + 1 == sample_steps[sample_idx]) {
            double E = T(p, dof) + V(q, dof);
            double Er = fabs(E - H_ref) / H_ref;
            lle[sample_idx] = sum_log / time;
            if (fp) {
                fprintf(fp, "%.16f %.16e %.16e %.16e\n", time, Er,
                        lle[sample_idx], standard_error);
                fflush(fp);
            }
            sample_idx++;
        }
    }

    free(X);
    free(dV);
    free(q2);
    free(Hdq);
    free(lle_vals);
}
