#ifndef LINALG_H_
#define LINALG_H_

//! \file linalg.h

#include <cmath>

//!Crout forward substitution
/*!
Performs forward substitution to solve a lower triangular matrix equation where the elements on the main diagonal of the lower triangular matrix are assumed to all be 1.
\param[in] L n by n lower triangular matrix
\param[in] b right hand side of matrix equation U*out = b
\param[in] p array of permutation indices from crout_LU
\param[in] n dimensions of A and length of b
\param[out] out array of length n, result of substitution
*/
void crout_forw_sub (double **L, double *b, int *p, int n, double *out);

//!Backward substitution
/*!
Performs backward substitution to solve an upper triangular matrix equation
input:
\param[in] U n by n upper triangular matrix
\param[in] b right hand side of matrix equation U*out = b
\param[in] n dimensions of A and length of b
\param[out] out array of length n, result of substitution
*/
void back_sub (double **U, double *b, int n, double *out);

//!Crout LU decomposition
/*!
performs LU decomposition with partial pivoting on the matrix A and sets permutation indices in the array p, following the outline in the book:
    + Demmel, James W. Applied numerical linear algebra. Vol. 56. Siam, 1997.
\param A matrix of size n to decompose (overwritten)
\param[in] n size of square matrix
\param[out] p integer array loaded with permutation indices
\return success code, which is 0 for success and 1 for failure due to singular matrix
*/
int crout_LU (double **A, int n, int *p);

//!Solves a matrix equation where the matrix has already be crout LU decomposed
/*!
\param[in] LU crout LU decomposed matrix
\param[in] p permutation indices
\param[in] b right-hand side of matrix equation
\param[in] n size of matrix and arrays
\param[out] out solution to linear system
*/
void solve_LU (double **LU, int *p, double *b, int n, double *out);

#endif
