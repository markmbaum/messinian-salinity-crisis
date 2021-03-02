//! \file msc_gcv.cc

#include "msc_gcv.h"

MscGcv::MscGcv () :
    OdeVern65 (2), //ode system of two equations
    Newton (2),    //nonlinear Newton system of the same two equations
    c1 (2.068e12), //fit parameter for fAm() and fzo()
    a1 (2754),     //fit parameter for fAm() and fzo()
    c2 (4.035e11), //fit parameter for fAm() and fzo()
    a2 (127.5),    //fit parameter for fAm() and fzo()
    Ao (360.0e12) {//area of world ocean without Mediterranean

    //default system name
    set_name("msc_gcv");

    //default initial conditions
    set_sol(0, -60.0); //sill level
    set_sol(0, 0.0); //Mediterranean level

    //default properties of built-in adaptive time step selection algorithm
    set_tol(1e-6);
    set_facmin(1e-2);
    set_facmax(1e1);
}

//------------------------------------------------------------------------------
//Mediterranean functions

double MscGcv::fAm (double zm) {
    return( c1*exp(zm/a1) + c2*exp(zm/a2) );
}

double MscGcv::fzo (double zm) {
    return( (c1*a1/Ao)*(1 - exp(zm/a1)) + (c2*a2/Ao)*(1 - exp(zm/a2)) );
}

//------------------------------------------------------------------------------
//system of ODEs

void MscGcv::ode_fun (double *solin, double *fout) {

    double zs, zm, zo, Am, tau, S, Q, T;
    //label incoming sill and Mediterranean levels
    zs = solin[0];
    zm = solin[1];
    //Mediterranean area and ocean level are functions of zm
    Am = fAm(zm);
    zo = fzo(zm);
    //channel slope
    /*NOTE:
    The slope is not allowed to be negative here for numerical reasons.
    If the model behaves physically, the slope should never be negative because
    the flow velocity over the sill goes to zero at exactly zo = zm. The basin
    should never fill above the ocean, not even by a tiny amount. In some cases,
    however, a tiny negative slope appears somewhere in the time stepping
    process. This is almost certainly just a numerical inaccuracy that occurs
    when the Mediterranean is filling extremely rapidly. However, if S is ever
    a negative number, even a vanishingly small one, it blows up the
    calculation of Q upon the evaluation of S^(13/14). When a NaN appears
    there, it carries through all subsequent calculations and the model fails.
    So, even though it should not be necessary on a physical basis, negative
    slope values are forced to zero for numerical feasibility. This issue
    cropped up in lots of different ODE solving algorithms, explicity and
    implicit, so it is not specific to the algorithm used here.*/
    S = max((zo - zm)/L, 0.0);
    //sill shear stress
    tau = RHO*GRAV*(zo - zs)*S;
    //bundle of constants
    T = Cw*pow((tauc + U/kb)/(RHO*GRAV), -3.0/13);
    //discharge for Turowski width
    Q = (pow(T, 13.0/7)/n)*pow(max(zo - zs, 0.0), 65.0/21)*pow(S, 13.0/14);
    //sill level time derivative, accounting for the critical τ threshold
    fout[0] = (tau > tauc) ? U - kb*pow(tau - tauc, a) : U;
    //Mediterranean level time derivative
    fout[1] = P - E + (R + Q)/Am;
}

//------------------------------------------------------------------------------

int MscGcv::has_root (double zslo, double zshi, double zmlo, double zmhi,
                      double *rzs, double *rzm, int nzs, int nzm) {

    //set the root values to NaN initially
    *rzs = NAN; *rzm = NAN;
    //success code
    int suc = -1;
    //get vectors of zs and zm values
    std::vector<double> zs = linspace(zslo, zshi, nzs);
    std::vector<double> zm = linspace(zmlo, zmhi, nzm);
    //allocate an array for the root
    double x[2];
    //sweep the grid for roots
    for (int i=0; i<nzs; i++) {
        for (int j=0; j<nzm; j++) {
            //initial values
            x[0] = zs[i]; x[1] = zm[j];
            //attempt to find a root
            suc = solve_Newton(x);
            //if successful, set the output and return a success code
            if ( suc == 0 ) {
                *rzs = x[0]; *rzm = x[1];
                return(suc);
            }
        }
    }
    //return the previous failure code
    return(suc);
}

int MscGcv::has_oscillation (double zs0, double zm0,
                             double tint, double tlim,
                             double rootrate) {

    //labels
    double zs, zm, dzs, dzm;
    double fout[2];
    double tprev, zsprev, zmprev;
    //root finding variables
    int rsuc;
    double r[2];

    //store current state
    tprev = get_t();
    zsprev = get_sol(0);
    zmprev = get_sol(1);

    //reset time and set initial conditions
    set_t(0.0);
    set_sol(0, zs0);
    set_sol(1, zm0);

    int osc = 0; //initialize to zero, which indicates oscillations
    //small initial solve to initialize adaptive time step size
    solve_adaptive(YRSEC, YRSEC/100, false);
    while ( get_t() < tlim ) {
        //integrate for a time interval of tint
        solve_adaptive(tint, get_dt(), false);
        //label the current solution
        zs = get_sol(0);
        zm = get_sol(1);
        //get the current derivatives
        ode_fun(get_sol(), fout);
        dzs = fout[0];
        dzm = fout[1];
        //if the sill is above the ocean, it will be forever
        if ( zs > fzo(zm) ) {
            osc = 1;
            break;
        }
        //quick fixed point check against 1 micron/year rate of change
        ode_fun(get_sol(), fout); //evaluate current derivatives
        if ( (fabs(dzs) < rootrate) && (fabs(dzm) < rootrate) ) {
            //try to find a root near the current state
            r[0] = zs; r[1] = zm;
            rsuc = solve_Newton(r);
            //check if the system is very near a fixed pt, if one was found
            if ( rsuc == 0 ) {
                if ( is_close(r[0], zs) && is_close(r[1], zm) ) {
                    osc = -1;
                    break;
                }
            }
        }

    }

    //put things back where they were found
    set_t(tprev);
    set_sol(0, zsprev);
    set_sol(1, zmprev);

    //return the oscillation code
    //  1 = sill above ocean, no oscillation
    //  0 = oscillation
    // -1 = fixed pt with sill below ocean, steady flow, no oscillation
    return(osc);
}

double MscGcv::index (double zs, double zm, double r) {

    double tol = 3.0*PI/4.0;
    //array for angle of vector field
    int n = 100000;
    double phi[100000];
    //calculate angles
    double theta;
    double z[2];
    double v[2];
    for (long i=0; i<n; i++) {
        //angle around circle surrounding x0,y0
        theta = 2*PI*double(i)/double(n-1);
        //coordinates of sample point
        z[0] = zs + r*cos(theta);
        z[1] = zm + r*sin(theta);
        //vector components
        ode_fun(z, v);
        //angle of components
        phi[i] = atan2(v[1], v[0]);
    }
    //correct for jumps
    double jump;
    for (long i=1; i<n; i++) {
        if ( (phi[i] > tol && phi[i-1] < -tol) || (phi[i] < -tol && phi[i-1] > tol) ) {
            jump = (phi[i-1] > 0.0) ? 2*PI : -2*PI;
            for (long j=i; j<n; j++) phi[j] += jump;
        }
    }

    return( (phi[n-1] - phi[0])/(2*PI) );
}

void MscGcv::print () {
    printf("zs = %-+g m\n", get_sol(0));
    printf("zm = %-+g m\n", get_sol(1));
    printf("zo = %-+g m\n", fzo(get_sol(1)));
}
