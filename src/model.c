#include "../include/model.h"
#include "../include/vectors.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double V(const double *q, const int dof) {
    double S1 = 0.0, S2 = 0.0, S3 = 0.0, S4 = 0.0;
    double sum_single = 0.0;

    // Compute sums and single-particle energy
    for (int i = 0; i < dof; i++) {
        double qi = q[i];
        double qi2 = qi * qi;
        double qi3 = qi2 * qi;
        double qi4 = qi3 * qi;
        S1 += qi;
        S2 += qi2;
        S3 += qi3;
        S4 += qi4;

        sum_single += 0.5 * qi2 + 0.25 * qi4;
    }

    // Interaction contribution (O(N))
    double sum_inter =
        (dof * S4 - 4 * S3 * S1 + 3 * S2 * S2) / (4.0 * (dof - 1));

    return sum_single + sum_inter;
}

void dVdq(double *dVdq_out, double const *q, const int dof) {

    double factor = 1.0 / (dof - 1);

    double S1 = 0.0;
    double S2 = 0.0;
    double S3 = 0.0;
    for (int i = 0; i < dof; i++) {
        double xi = q[i];
        S1 += xi;
        S2 += xi * xi;
        S3 += xi * xi * xi;
    }

    for (int n = 0; n < dof; n++) {
        double qn = q[n];
        double qn2 = qn * qn;
        double qn3 = qn2 * qn;
        double coup = dof * qn3 - 3.0 * qn2 * S1 + 3.0 * qn * S2 - S3;
        dVdq_out[n] = qn + qn3 + factor * coup;
    }
}

void hessian_V(double *H, const double *q, double *q2, const int dof) {

    // Precompute S1 = Σ q, S2 = Σ q²
    double S1 = 0.0;
    double S2 = 0.0;

    for (int i = 0; i < dof; i++) {
        double qi = q[i];
        S1 += qi;
        S2 += qi * qi;
        q2[i] = qi * qi;
    }

    // Fill the matrix column major
    for (int i = 0; i < dof; i++) {
        double qi = q[i];
        double q2i = q2[i];
        H[i + i * dof] = 1.0 + 3.0 * q2i +
                         3.0 * (dof * q2i - 2.0 * qi * S1 + S2) / (dof - 1.0);
        for (int j = 0; j < dof; j++) {
            if (i != j) {
                double qj = q[j];
                double q2j = q2[j];
                double diff2 = q2i - 2.0 * qi * qj + q2j;
                H[i + j * dof] = -3.0 * diff2 / (dof - 1.0);
            }
        }
    }
}

void hessian_V_times_vec(double *y, const double *q, const double *q2,
                         size_t dof, double S1, double S2, const double *v) {
    const double denom = (double)dof - 1.0;
    const double inv = 1.0 / denom;

    // Sv, Sqv, Sq2v
    double Sv = 0.0, Sqv = 0.0, Sq2v = 0.0;
    for (size_t i = 0; i < dof; i++) {
        const double vi = v[i];
        Sv += vi;
        Sqv += q[i] * vi;
        Sq2v += q2[i] * vi;
    }

    // y_n formula
    for (size_t n = 0; n < dof; n++) {
        const double qn = q[n];
        const double qn2 = q2[n];
        const double vn = v[n];

        const double coeff_vn =
            1.0 + 3.0 * qn2 +
            inv * (3.0 * (double)dof * qn2 - 6.0 * qn * S1 + 3.0 * S2);

        const double meanfield =
            (3.0 * inv) * (qn2 * Sv - 2.0 * qn * Sqv + Sq2v);

        y[n] = vn * coeff_vn - meanfield;
    }
}

double T(double const *p, const int dof) {
    double sum_p2 = 0.0;
    for (int i = 0; i < dof; i++)
        sum_p2 += p[i] * p[i];
    return sum_p2 / 2.0;
}

void calculate_eff_site_energies(const double *q, const double *p, size_t dof,
                                 double *eff_site_energy) {

    double S[4] = {0.0, 0.0, 0.0, 0.0};
    double mu[3] = {0.0, 0.0, 0.0};
    for (size_t i = 0; i < dof; i++) {
        double qi = q[i];
        double qi2 = qi * qi;
        double qi3 = qi2 * qi;
        double qi4 = qi3 * qi;

        S[0] += qi;
        S[1] += qi2;
        S[2] += qi3;
        S[3] += qi4;
    }

    for (size_t i = 0; i < 3; i++)
        mu[i] = S[i] / dof;

    for (size_t i = 0; i < dof; i++) {
        double qi = q[i];
        double qi2 = qi * qi;
        double qi3 = qi2 * qi;
        double qi4 = qi3 * qi;
        double pi = p[i];
        double pi2 = pi * pi;

        double V_eff = 0.5 * qi2 + 0.5 * qi4 - mu[0] * qi3 +
                       3 * mu[1] * qi2 / 2 - mu[2] * qi;

        eff_site_energy[i] = 0.5 * pi2 + V_eff;
    }
}

void generate_waterbag_initial_condition(size_t dof, double eps, double *q,
                                         double *p, unsigned int seed) {

    for (size_t i = 0; i < dof; i++) {
        q[i] = 0;
        double rnd = (double)rand_r(&seed) / (double)RAND_MAX;
        p[i] = rnd - 0.5;
    }

    double sum_p = sum_vector(p, dof);
    for (size_t i = 0; i < dof; i++)
        p[i] -= sum_p / dof;

    double total_energy = T(p, dof);
    double target_energy = eps * dof;

    if (total_energy <= 0) {
        fprintf(
            stderr,
            "Error: nonpositive total energy in generate_initial_condition\n");
        exit(EXIT_FAILURE);
    }

    double scale = sqrt(target_energy / total_energy);
    for (size_t i = 0; i < dof; i++)
        p[i] *= scale;

    total_energy = T(p, dof) + V(q, dof);
    double diff = fabs(total_energy - target_energy);
    double rel_error = diff / target_energy;
    if (rel_error > 1e-10) {
        fprintf(stderr, "Error: target_energy not reached\n");
        fprintf(stderr, "total_energy = %.16f\ntarget_energy = %.16f\n",
                total_energy, target_energy);
        exit(EXIT_FAILURE);
    }
}
