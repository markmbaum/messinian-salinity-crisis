//! \file linalg.cc

#include "linalg.h"

void crout_forw_sub (double **L, double *b, int *p, int n, double *out) {

    for (int i=0; i<n; i++) {
        out[i] = b[p[i]]; //the permutation is handled here
        for (int j=0; j<i; j++) {
            out[i] -= L[i][j]*out[j];
        }
        //element on the diagonal is 1 for crout LU decomposed matrices
    }
}

void back_sub (double **U, double *b, int n, double *out) {

    for (int i=n-1; i>=0; i--) {
        out[i] = b[i];
        for (int j=i+1; j<n; j++) {
            out[i] -= U[i][j]*out[j];
        }
        out[i] /= U[i][i];
    }
}

int crout_LU (double **A, int n, int *p) {

    int i, j, k, idx, ti;
    double m, td;

    //initialize the permutation array
    for (i=0; i<n; i++) p[i] = i;

    //move through the matrix
    for (i=0; i<n; i++) {
        //find the largest element in the column beneath A[i][i]
        m = 0.0;
        idx = i;
        for (j=i; j<n; j++) {
            td = fabs(A[j][i]);
            if ( td > m ) {
                m = td;
                idx = j;
            }
        }
        //if all elements are zero, the matrix is singular
        if ( !(fabs(m) > 0.0) ) return(1); //failure!
        //if the largest element isn't on the diagonal, pivot
        if ( idx != i ) {
            //swap rows
            for (j=0; j<n; j++) {
                td = A[i][j];
                A[i][j] = A[idx][j];
                A[idx][j] = td;
            }
            //swap permutation indices
            ti = p[i];
            p[i] = p[idx];
            p[idx] = ti;
        }
        //continue with the decomposition
        for (j=i+1; j<n; j++)
            A[j][i] /= A[i][i];
        for (j=i+1; j<n; j++)
            for (k=i+1; k<n; k++)
                A[j][k] -= A[j][i]*A[i][k];
    }
    //success!
    return(0);
}

void solve_LU (double **LU, int *p, double *b, int n, double *out) {

    //run forward substitution, handling the permutation along the way
    crout_forw_sub(LU, b, p, n, out);
    //run backward substitution, overwriting "out" along the way
    back_sub(LU, out, n, out); //solution is now in "out"
}
