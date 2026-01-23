#include "../include/integrator.h"
#include "../include/model.h"
#include <stddef.h>
#include <stdlib.h>

static const double ALPHA =
    1.0 / (2.0 - 1.25992104989487); // ~0.6756035959798289
static const double BETA =
    -1.25992104989487 / (2.0 - 1.25992104989487); // ~-0.3512071919596578

void symplectic_step(double *q, double *p, const size_t dof, double time_step,
                     double *dV) {

    // Half kick
    dVdq(dV, q, dof);
    for (size_t i = 0; i < dof; i++)
        p[i] -= 0.5 * time_step * dV[i];

    // Drift
    for (size_t i = 0; i < dof; i++)
        q[i] += time_step * p[i];

    // Half kick
    dVdq(dV, q, dof);
    for (size_t i = 0; i < dof; i++)
        p[i] -= 0.5 * time_step * dV[i];
}

void symplectic_step_traj_tan(double *q, double *p, double *q2, double *X,
                              double *dV, double *Hdq, size_t dof, size_t n_dev,
                              double time_step) {

    // --- Half kick on main trajectory ---
    dVdq(dV, q, dof); // gradient at current q
    for (size_t i = 0; i < dof; i++)
        p[i] -= 0.5 * time_step * dV[i];

    // --- Half kick on tangent momenta ---
    double S1, S2;
    precompute_S1_S2_q2(q, q2, dof, &S1, &S2);
    for (size_t k = 0; k < n_dev; k++) {
        double *Xk = X + k * 2 * dof;
        double *dq = Xk;
        double *dp = Xk + dof;

        hessian_V_times_vec(Hdq, q, q2, dof, S1, S2, dq);

        for (size_t i = 0; i < dof; i++)
            dp[i] -= 0.5 * time_step * Hdq[i];
    }

    // --- Drift ---
    for (size_t i = 0; i < dof; i++)
        q[i] += time_step * p[i];

    // --- Drift on tangent positions ---
    for (size_t k = 0; k < n_dev; k++) {
        double *Xk = X + k * 2 * dof;
        double *dq = Xk;
        double *dp = Xk + dof;
        for (size_t i = 0; i < dof; i++)
            dq[i] += time_step * dp[i];
    }

    // --- Half kick on main trajectory ---
    dVdq(dV, q, dof);
    for (size_t i = 0; i < dof; i++)
        p[i] -= 0.5 * time_step * dV[i];

    // --- Half kick on tangent momenta ---
    precompute_S1_S2_q2(q, q2, dof, &S1, &S2);
    for (size_t k = 0; k < n_dev; k++) {
        double *Xk = X + k * 2 * dof;
        double *dq = Xk;
        double *dp = Xk + dof;

        hessian_V_times_vec(Hdq, q, q2, dof, S1, S2, dq);

        for (size_t i = 0; i < dof; i++)
            dp[i] -= 0.5 * time_step * Hdq[i];
    }
}

void yoshida4_step(double *q, double *p, const size_t dof, double *dV) {
    symplectic_step(q, p, dof, ALPHA * TIME_STEP, dV);
    symplectic_step(q, p, dof, BETA * TIME_STEP, dV);
    symplectic_step(q, p, dof, ALPHA * TIME_STEP, dV);
}

void yoshida4_step_traj_tan(double *q, double *p, double *q2, double *X,
                            double *dV, double *Hdq, size_t dof, size_t n_dev) {
    // ALPHA step
    symplectic_step_traj_tan(q, p, q2, X, dV, Hdq, dof, n_dev,
                             ALPHA * TIME_STEP);
    // BETA step
    symplectic_step_traj_tan(q, p, q2, X, dV, Hdq, dof, n_dev,
                             BETA * TIME_STEP);
    // ALPHA step
    symplectic_step_traj_tan(q, p, q2, X, dV, Hdq, dof, n_dev,
                             ALPHA * TIME_STEP);
}
