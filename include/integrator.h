#pragma once

#include <stddef.h>

#define TIME_STEP 0.01

/* Symplectic step (trajectory) */
void symplectic_step(double *q, double *p, const size_t dof, double time_step,
                     double *dV);

/* Symplectic step (trajectory + tangent) */
void symplectic_step_traj_tan(double *q, double *p, double *q2, double *X,
                              double *dV, double *Hv, size_t dof, size_t n_dev,
                              double time_step);

/* 4th-order Yoshida integration (trajectory) */
void yoshida4_step(double *q, double *p, const size_t dof, double *dV);

/* 4th-order Yoshida integration (trajectory + tangent) */
void yoshida4_step_traj_tan(double *q, double *p, double *q2, double *X,
                            double *dV, double *Hv, size_t dof, size_t n_dev);
