#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/allocs.h"

// Fortran LAPACK routines
extern void dgeqrf_(const int *m, const int *n, double *a, const int *lda,
                    double *tau, double *work, const int *lwork, int *info);

extern void dorgqr_(const int *m, const int *n, const int *k, double *a,
                    const int *lda, double *tau, double *work, const int *lwork,
                    int *info);

// Fortran LAPACK routine for SVD
extern void dgesvd_(const char *JOBU, const char *JOBVT, const int *M,
                    const int *N, double *A, const int *LDA, double *S,
                    double *U, const int *LDU, double *VT, const int *LDVT,
                    double *WORK, const int *LWORK, int *INFO);

// Transpose helpers
void transpose_row_col(const double *in, double *out, int m, int n) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            out[j * m + i] = in[i * n + j]; // row-major -> column-major
}

void transpose_col_row(const double *in, double *out, int m, int n) {
    for (int j = 0; j < n; j++)
        for (int i = 0; i < m; i++)
            out[i * n + j] = in[j * m + i]; // column-major -> row-major
}

// QR factorization for column-major input/output
// Input: A_col (m x n) column-major
// Output: Q in A_col (m x n) column-major, R (n x n) column-major in R_col
int qr_colmajor(int m, int n, double *A_col, double *R_col) {
    int info;
    int lda = m; // leading dimension for column-major

    // Allocate tau
    double *tau = malloc(n * sizeof(*tau));
    if (!tau) {
        perror("tau allocation failed");
        return -1;
    }

    // Workspace query
    int lwork = -1;
    double wkopt;
    dgeqrf_(&m, &n, A_col, &lda, tau, &wkopt, &lwork, &info);
    if (info != 0) {
        free(tau);
        return info;
    }
    lwork = (int)wkopt;
    double *work = malloc(lwork * sizeof(*work));
    if (!work) {
        perror("work allocation failed");
        free(tau);
        return -1;
    }

    // Compute QR factorization
    dgeqrf_(&m, &n, A_col, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        free(work);
        free(tau);
        return info;
    }

    // Extract R from A_col (column-major)
    // R is upper-triangular n x n
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double val = (i <= j) ? A_col[j * lda + i] : 0.0; // column-major
            R_col[i + j * n] = val; // store R in column-major
        }
    }

    // Generate Q explicitly in column-major
    dorgqr_(&m, &n, &n, A_col, &lda, tau, work, &lwork, &info);

    free(work);
    free(tau);

    return info;
}

// QR factorization with row-major input/output
// Input: A (m x n) row-major
// Output: Q (m x n) row-major in A, R (n x n) row-major in R
int qr_rowmajor(int m, int n, double *A_row, double *R_row) {
    int info;
    int lda = m; // leading dimension for column-major

    // Allocate column-major copy of A
    double *A_col = malloc(m * n * sizeof(*A_col));
    if (!A_col) {
        perror("A_col allocation failed");
        return -1;
    }
    transpose_row_col(A_row, A_col, m, n);

    // Allocate tau
    double *tau = malloc(n * sizeof(*tau));
    if (!tau) {
        perror("tau allocation failed");
        free(A_col);
        return -1;
    }

    // Workspace query
    int lwork = -1;
    double wkopt;
    dgeqrf_(&m, &n, A_col, &lda, tau, &wkopt, &lwork, &info);
    if (info != 0) {
        free(A_col);
        free(tau);
        return info;
    }
    lwork = (int)wkopt;
    double *work = malloc(lwork * sizeof(*work));
    if (!work) {
        perror("work allocation failed");
        free(A_col);
        free(tau);
        return -1;
    }

    // Compute QR factorization
    dgeqrf_(&m, &n, A_col, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        free(work);
        free(tau);
        free(A_col);
        return info;
    }

    // Extract R from column-major A_col and transpose to row-major R_row
    // R is n x n upper triangle
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            double val =
                (i <= j) ? A_col[j * lda + i] : 0.0; // A(i,j) col-major
            R_row[i * n + j] = val; // directly row-major in R_row
        }

    // Generate Q explicitly in column-major
    dorgqr_(&m, &n, &n, A_col, &lda, tau, work, &lwork, &info);
    if (info != 0) {
        free(work);
        free(tau);
        free(A_col);
        return info;
    }

    // Transpose Q back to row-major into A_row
    transpose_col_row(A_col, A_row, m, n);

    free(work);
    free(tau);
    free(A_col);

    return info;
}

// Compute singular values of a column-major matrix A (m x n)
// Input:  A_col: column-major m x n
//         m, n: dimensions
// Output: S: array of min(m,n) singular values (descending order)
int singular_values_colmajor(int m, int n, double *A, double *S) {
    int info;
    int lda = m;
    int ldu = 1;  // since JOBU = 'N', LDU can be 1
    int ldvt = 1; // since JOBVT = 'N', LDVT can be 1

    char jobu = 'N';  // don't compute U
    char jobvt = 'N'; // don't compute V^T

    double *A_col = xmalloc(m * n, sizeof *A_col);
    memcpy(A_col, A, m * n * sizeof(*A_col));

    // Workspace query
    int lwork = -1;
    double wkopt;
    dgesvd_(&jobu, &jobvt, &m, &n, A_col, &lda, S, NULL, &ldu, NULL, &ldvt,
            &wkopt, &lwork, &info);
    if (info != 0)
        return info;

    lwork = (int)wkopt;
    double *work = xmalloc(lwork, sizeof *work);
    if (!work) {
        perror("work allocation failed");
        return -1;
    }

    // Compute singular values
    dgesvd_(&jobu, &jobvt, &m, &n, A_col, &lda, S, NULL, &ldu, NULL, &ldvt,
            work, &lwork, &info);

    free(work);
    free(A_col);
    return info;
}
