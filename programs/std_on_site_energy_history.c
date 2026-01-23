#include "../include/allocs.h"
#include "../include/integrator.h"
#include "../include/math_helpers.h"
#include "../include/model.h"
#include "../include/parse_args.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define TIME_WINDOW 500.0

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
    size_t window_steps =
        (size_t)(TIME_WINDOW / TIME_STEP); // Number of window steps

    double *dV =
        xmalloc(dof, sizeof *dV); // Preallocated array for the gradient of V
    double *q = xmalloc(dof, sizeof *q); // Coordinates
    double *p = xmalloc(dof, sizeof *p); // Momenta
    generate_waterbag_initial_condition(dof, eps, q, p, seed);

    char filename[512];
    // Build the filename from the parameters
    snprintf(
        filename, sizeof(filename),
        "data/std_on_site_energy_history_dof=%zu_eps=%.6f_T=%.1f_seed=%d.dat",
        dof, eps, total_time, seed);
    // Open the file
    FILE *fp = fopen(filename, "w");
    if (!fp) {
        perror("fopen");
        exit(EXIT_FAILURE);
    }

    // Reference energy to check numerical integration
    double H_ref = T(p, dof) + V(q, dof);

    // Array for the effective site energy
    double *all_eff_site_energy = xmalloc(dof, sizeof *all_eff_site_energy);

    // Arrays for the Welford's standard deviation calculation
    double *mean = xcalloc(dof, sizeof *mean);
    double *M2 = xcalloc(dof, sizeof *M2);

    size_t count = 0;
    for (size_t n = 0; n < n_steps; n++) {

        yoshida4_step(q, p, dof, dV);

        calculate_eff_site_energies(q, p, dof, all_eff_site_energy);

        welford_update(all_eff_site_energy, mean, M2, dof, &count);

        if ((n + 1) % window_steps == 0) {

            double H = T(p, dof) + V(q, dof);
            double Er = fabs(H - H_ref) / H_ref;

            double sigma = welford_std(M2, dof, count);

            fprintf(fp, "%.16f %.16e %.16e\n", (n + 1) * TIME_STEP, Er, sigma);
            fflush(fp);

            // Reset for the new window
            count = 0;
            welford_reset(mean, M2, dof);
        }
    }

    free(dV);
    free(q);
    free(p);
    free(all_eff_site_energy);
    free(mean);
    free(M2);
    fclose(fp);
    return 0;
}
