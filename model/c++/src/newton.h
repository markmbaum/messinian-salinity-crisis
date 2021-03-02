#ifndef NEWTON_H_
#define NEWTON_H_

//! \file newton.h

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include "linalg.h"

//!Newton's method for nonlinear systems of equations
/*!
This class implements Newton's method for nonlinear systems of equations. The virtual functions F_Newton and Jac_Newton allow a derived class to implement the system of equations and its Jacobian matrix.
*/
class Newton {

    public:
        //!constructs
        /*!
        \param[in] n size of system
        */
        Newton (unsigned long n);
        //!destructs
        virtual ~Newton ();

        //-------------------
        //getters and setters

        //!gets the size of the system
        unsigned long get_n () { return(n_); }
        //!gets the L infinity tolerance
        double get_tol_Newton () { return(tol_Newton_); }
        //!gets iteration counter
        unsigned long get_iter_Newton () { return(iter_Newton_); }
        //!gets the LU decomposition interval
        int get_iJLU () { return(iJLU_); }
        //!gets the LU decomposition counter
        unsigned long get_nJLU () { return(nJLU_); }
        //!gets the LU solve counter
        unsigned long get_n_solve_LU () { return(n_solve_LU_); }
        //!gets whether modified Newtion's is being used
        bool get_modified () { return(modified_); }
        //!gets whether no LU decompositions should be done
        bool get_ignore_JLU () { return(ignore_JLU_); }
        //!gets adjustment factor for numerical jacobian
        double get_absjacdel () { return(absjacdel_); }
        //!gets adjustment factor for numerical jacobian
        double get_reljacdel () { return(reljacdel_); }

        //!sets the L infinity tolerance
        void set_tol_Newton (double tol_Newton) { tol_Newton_ = tol_Newton; }
        //!sets the  iteration counter
        void set_iter_Newton (unsigned long iter_Newton) { iter_Newton_ = iter_Newton; }
        //!sets the LU decomposition interval
        void set_iJLU (int iJLU) { iJLU_ = iJLU; }
        //!sets whether modified Newtion's is being used
        void set_modified (bool modified) { modified_ = modified; }
        //!sets whether no LU decompositions should be done
        void set_ignore_JLU (bool ignore_JLU) { ignore_JLU_ = ignore_JLU; }
        //!sets adjustment factor for numerical jacobian
        void set_absjacdel (double absjacdel) { absjacdel_ = absjacdel; }
        //!sets adjustment factor for numerical jacobian
        void set_reljacdel (double rbsjacdel) { reljacdel_ = rbsjacdel; }

        //!Solve the system of equations
        /*!Solve for a root of the function f_Newton using the jacobian J_Newton, both of which must be implemented in derived classes
        \param x initial guess for the root and final value at root
        \return success code, which takes values
            0. success!
            1. failure...singular matrix encountered
            2. failure...too many iterations
            3. failure...solution contains nan(s)
            4. failure...solution contains inf(s)*/
        int solve_Newton (double *x);

    protected:

        //!evaluates the function being zeroed
        /*!
        \param[in] x current values of system
        \param[out] f evaluated system of equation
        */
        virtual void f_Newton (double *x, double *f) = 0;

        //!evaluates the Jacobian matrix of the function being zeroed
        /*!
        Unless implemented by the derived class the Jacobian will be estimated with a hasty finite differences algorithm
        \param[in] x current values of system
        \param[out] J evaluated Jacobian of the system
        */
        virtual void J_Newton (double *x, double **J);

    private:

        //size of system
        unsigned long n_;
        //whether to use only a single evaluation and LU decomposition of Jac
        bool modified_;
        //whether to do no JLU updates at all
        bool ignore_JLU_;
        //error tolerance
        double tol_Newton_;
        //iteration tolerance
        unsigned long iter_Newton_;
        //number of iterations after which the Jacobian is updated and redecomposed
        int iJLU_;
        //number of times Jac has been updated
        unsigned long nJLU_;
        //number of times the LU decomposed matrix has been used to solve a matrix eq
        unsigned long n_solve_LU_;
        //iteration interval for checking solution integrity
        unsigned long icheck_;
        //current evaluation of F
        double *f_;
        //current evaluation of J
        double **J_;
        //update array
        double *delx_;
        //permutation array
        int *p_;
        //!absolute adjustment fraction for numerical Jacobian, if needed
        double absjacdel_;
        //!relative adjustment fraction for numerical Jacobian, if needed
        double reljacdel_;
        //an extra array for evaluating numerical jacobian
        double *g_;

        //finds the infinity norm of the update vector delx and of y
        void err (double *errx, double *erry);
        //recomputes the Jacobian and crout LU decomposes it
        int JLU (double *x);
        //wrapper around LU solving routine for counting solves
        void solve_LU_(double **LU, int *p, double *b, int n, double *out);
        //function for checking solution integrity (can't have NAN or INFINITY)
        int check_integrity (double *x);
};

#endif
