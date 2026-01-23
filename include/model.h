#pragma once

#include <stddef.h>

/* Potential energy */
double V(double const *q, const int dof);

/* Gradient of potential */
void dVdq(double *dVdq_out, double const *q, const int dof);

/* Hessian of potential */
void hessian_V(double *H, const double *q, double *q2, const int dof);

void hessian_V_times_vec(double *y, const double *q, const double *q2,
                         size_t dof, double S1, double S2, const double *v);

/* Kinetic energy */
double T(double const *p, const int dof);

static inline void precompute_S1_S2_q2(const double *q, double *q2, size_t dof,
                                       double *S1, double *S2) {
    double s1 = 0.0, s2 = 0.0;
    for (size_t i = 0; i < dof; i++) {
        const double qi = q[i];
        const double qi2 = qi * qi;
        q2[i] = qi2;
        s1 += qi;
        s2 += qi2;
    }
    *S1 = s1;
    *S2 = s2;
}

void calculate_eff_site_energies(const double *q, const double *p, size_t dof,
                                 double *eff_site_energy);

void generate_waterbag_initial_condition(size_t dof, double eps, double *q,
                                         double *p, unsigned int seed);
