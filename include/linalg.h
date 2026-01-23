#ifndef LINALG_H
#define LINALG_H

// QR factorization (row-major interface)
// Input: A (m x n) row-major
// Output: Q in A (row-major), R (n x n) row-major
// Returns: 0 on success, LAPACK error code otherwise
int qr_rowmajor(int m, int n, double *A, double *R);

int qr_colmajor(int m, int n, double *A_col, double *R_col);

int singular_values_colmajor(int m, int n, double *A_col, double *S);

#endif // LINALG_H
