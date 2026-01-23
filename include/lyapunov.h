#pragma once

#include <stddef.h>
#include <stdio.h>

void LLE_history_stop_condition(double *q, double *p, size_t dof,
                                double total_time, size_t *sample_steps,
                                double *lle, double rtol, double atol,
                                size_t num_vals, FILE *fp);
