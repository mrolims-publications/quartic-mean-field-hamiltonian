#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/lyapunov.h"
#include "../include/math_helpers.h"
#include "../include/model.h"
#include "../include/parse_args.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char *argv[]) {

    if (argc != 8) {
        fprintf(
            stderr,
            "Usage: %s <dof> <eps> <total_time> <rtol> <atol> <num_vals> "
            "<seed>\n\n"
            "Arguments:\n"
            "  <dof>        Number of degrees of freedom (positive integer)\n"
            "  <eps>        Specific energy ε = E / N\n"
            "  <total_time> Maximum integration time\n"
            "  <rtol>       Relative tolerance for LLE convergence\n"
            "               (stop if std / mean < rtol over the last "
            "<num_vals>)\n"
            "  <atol>       Absolute tolerance for LLE convergence\n"
            "               (stop if std < atol over the last <num_vals>)\n"
            "  <num_vals>   Number of recent LLE values used in the "
            "convergence test\n"
            "  <seed>       Random number generator seed\n",
            argv[0]);
        return EXIT_FAILURE;
    }

    size_t dof = get_size_t(argc, argv, 1);
    double eps = get_double(argc, argv, 2);
    double total_time = get_double(argc, argv, 3);
    double rtol = get_double(argc, argv, 4);
    double atol = get_double(argc, argv, 5);
    size_t num_vals = get_size_t(argc, argv, 6);
    unsigned int seed = get_uint(argc, argv, 7);

    size_t n_steps = (size_t)(total_time / TIME_STEP);

    double *q = xmalloc(dof, sizeof *q);
    double *p = xmalloc(dof, sizeof *p);

    generate_waterbag_initial_condition(dof, eps, q, p, seed);

    char filename[512];
    snprintf(filename, sizeof(filename),
             "data/"
             "lle_history_dof=%zu_eps=%.6f_T=%.1f_dt=%.3f_rtol=%.16f_atol="
             "%.16f_num_vals=%zu_seed=%d.dat",
             dof, eps, total_time, TIME_STEP, rtol, atol, num_vals, seed);
    FILE *fp = fopen(filename, "w");

    size_t n_unique;
    size_t *sample_steps =
        ilogspace_unique(0.0, log10(n_steps), 10000, &n_unique);

    double *lle = xmalloc(n_unique, sizeof *lle);
    LLE_history_stop_condition(q, p, dof, total_time, sample_steps, lle, rtol,
                               atol, num_vals, fp);

    fclose(fp);
    free(q);
    free(p);
    free(sample_steps);
    free(lle);

    return 0;
}
