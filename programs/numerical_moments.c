#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/model.h"
#include "../include/parse_args.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define SAVE_EVERY 0.1

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

    size_t n_steps =
        (size_t)(total_time / TIME_STEP); // Number of integration steps
    size_t save_stride = (size_t)(SAVE_EVERY / TIME_STEP);

    double *dV =
        xmalloc(dof, sizeof *dV); // Preallocated array for the gradient of V
    double *q = xmalloc(dof, sizeof *q); // Coordinates
    double *p = xmalloc(dof, sizeof *p); // Momenta
    generate_waterbag_initial_condition(dof, eps, q, p, seed);

    char filename[512];
    // Build the filename from the parameters
    snprintf(filename, sizeof(filename),
             "data/moments_history_dof=%zu_eps=%.6f_T=%.1f_dt=%.4f_seed=%d.dat",
             dof, eps, total_time, TIME_STEP, seed);
    // Open the file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    for (size_t n = 0; n < n_steps; n++) {

        yoshida4_step(q, p, dof, dV);

        if ((n + 1) % save_stride == 0) {

            double time = (n + 1) * TIME_STEP;
            double M1 = 0.0;
            double M2 = 0.0;
            double M3 = 0.0;
            for (size_t i = 0; i < dof; i++) {
                double qi = q[i];
                double qi2 = qi * qi;
                double qi3 = qi2 * qi;
                M1 += qi;
                M2 += qi2;
                M3 += qi3;
            }

            M1 /= (double)dof;
            M2 /= (double)dof;
            M3 /= (double)dof;

            fprintf(fp, "%.16f %.16e %.16e %.16e\n", time, M1, M2, M3);
            fflush(fp);
        }
    }

    free(dV);
    free(q);
    free(p);
    fclose(fp);
    return 0;
}
