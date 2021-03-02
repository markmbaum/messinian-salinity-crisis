//! \file newton.cc

#include "newton.h"

Newton::Newton (unsigned long n) {

    //size of system
    n_ = n;
    //allocate an array for evaluating f
    f_ = new double[n];
    //allocate 2D array for evaluating J
    J_ = new double*[n];
    for (unsigned long i=0; i<n; i++) J_[i] = new double[n];
    //update array
    delx_ = new double[n];
    //permutation array
    p_ = new int[n];
    //error tolerance
    tol_Newton_ = 1e-8;
    //iteration limit
    iter_Newton_ = 100;
    //Jacobian update interval
    iJLU_ = 1;
    //JLU counter
    nJLU_ = 0;
    //number of times the LU decomposed matrix is used to solve a matrix eq
    n_solve_LU_ = 0;
    //iteration interval for checking solution integrity
    icheck_ = 2;
    //whether to use only a single evaluation and LU decomposition of Jac
    modified_ = false;
    //whether to do no Jac evaluations and LU decomposition at all
    ignore_JLU_ = false;
    //an extra array for evaluating numerical jacobian, if needed
    g_ = new double[n];
    //default adjustment factors for numerical jacobian
    absjacdel_ = 1e-7;
    reljacdel_ = 1e-7;
}

Newton::~Newton () {
    delete [] f_;
    for (unsigned long i=0; i<n_; i++) delete [] J_[i];
    delete [] J_;
    delete [] delx_;
    delete [] p_;
    delete [] g_;
}

void Newton::J_Newton (double *x, double **J) {

    unsigned long i,j;
    double delsol;

    //get derivatives at current position
    f_Newton(x, f_);
    //perturb and calculate finite differences
    for (i=0; i<n_; i++) {
        //perturb
        delsol = (absjacdel_ > x[i]*reljacdel_) ? absjacdel_ : x[i]*reljacdel_;
        x[i] += delsol;
        //evaluate function
        f_Newton(x, g_);
        //compute approximate derivatives
        for (j=0; j<n_; j++) J[j][i] = (g_[j] - f_[j])/delsol;
        //unperturb
        x[i] -= delsol;
    }
}

void Newton::err (double *errx, double *errf) {

    double tx, tf, mx = 0, my = 0;

    for (unsigned long i=0; i<n_; i++) {
        if ( (tx = fabs(delx_[i])) > mx ) mx = tx;
        if ( (tf = fabs(f_[i])) > my ) my = tf;
    }
    *errx = mx;
    *errf = my;
}

int Newton::JLU (double *x) {

    //evaluate the Jacobian matrix
    J_Newton(x, J_);
    //crout LU decompose the Jacobian, storing the result
    int suc = crout_LU(J_, n_, p_);
    //increment counter
    nJLU_++;
    //return the success code
    return(suc);
}

void Newton::solve_LU_ (double **LU, int *p, double *b, int n, double *out) {
    solve_LU(LU, p, b, n, out);
    n_solve_LU_++;
}

int Newton::check_integrity (double *x) {

    for (unsigned long i=0; i<n_; i++) {
        //any nans?
        if ( std::isnan(x[i]) ) return(3);
        //any infs?
        if ( !std::isfinite(x[i]) ) return(4);
    }
    return(0);
}

int Newton::solve_Newton (double *x) {

    unsigned long i, iter;
    int suc; //success (or failure) code, which is zero for success
    double errx, errf;

    //initialize the update vector
    for (i=0; i<n_; i++) delx_[i] = INFINITY;
    //initialize the error estimates
    err(&errx, &errf);

    //if modified Newton's, use only a single evaluation and LU decomposition of J
    if ( modified_ && (!ignore_JLU_) ) JLU(x);

    //iterate
    iter = 0;
    while ( (errx > tol_Newton_) || (errf > tol_Newton_) ) {
        //evaluate and LU decompose the Jacobian
        if ( (!ignore_JLU_) && (!modified_) && (iter % iJLU_ == 0) ) {
            suc = JLU(x);
            if (suc == 1) return(suc); //failure! singular matrix!
        }
        //evaluate the function
        f_Newton(x, f_);
        //swap the sign of f
        for (i=0; i<n_; i++) f_[i] = -f_[i];
        //solve the matrix equation
        solve_LU_(J_, p_, f_, n_, delx_);
        //update x
        for (i=0; i<n_; i++) x[i] += delx_[i];
        //update the error estimate
        err(&errx, &errf);
        //increment the counter
        iter++;
        //check iteration limits
        if ( iter > iter_Newton_ ) return(2);
        //check for nans and infs
        if ( iter % icheck_ == 0 ) {
            suc = check_integrity(x);
            if ( suc != 0 ) return(suc);
        }
    }

    return(0);
}
