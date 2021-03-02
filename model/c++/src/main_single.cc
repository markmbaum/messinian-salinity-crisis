//! \file main_single.cc

#include <string>

#include "omp.h"

#include "msc_gcv.h"


//!driver for running a single model integration
/*!
+ input parameters are set inside the function
+ output is written to the `out` directory
+ the output can be plotted with `scripts/plot_out.py`
*/
int main () {

    //output directory
    std::string dirout = "out";
    //integration time
    double tint = 1e2*KYRSEC;
    //initial time step
    double dt0 = YRSEC;

    //construct a model object
    MscGcv sys;

    //--------------
    //set parameters

    //erodability [m/yr/Pa^-a]
    sys.kb = 8e-6/YRSEC;
    //critical shear stress [Pa]
    sys.tauc = 50;
    //constant for Turowski channel width formula [?]
    sys.Cw = 6; //value not clear, but appears to be 6 in asalted.c/DEMO.sh??
    //uplift rate [m/s]
    sys.U = 4.9/MMYR;
    //erosion exponent []
    sys.a = 1.5;
    //stream length [m]
    sys.L = 100e3;

    //roughness coefficient [s/m^(1/3)]
    sys.n = 0.05;
    //precipitation rate [m/s]
    sys.P = 0.6/YRSEC;
    //evaporation rate [m/s]
    sys.E = 1.2/YRSEC;
    //input from various rivers [m^3/s]
    sys.R = 4500.0 + 12000.0;

    printf("\nsystem initialized\n");

    //------------------------
    //integrate system of ODEs

    //initial conditions
    sys.set_sol(0, -60); //sill level
    sys.set_sol(1, 0.0); //Mediterranean level

    printf("\nbeginning integration...\n");
    double tstart = omp_get_wtime();
    sys.solve_adaptive(tint, dt0, dirout.c_str(), 1);
    printf("integration finished after\n  %g seconds\n", omp_get_wtime() - tstart);
    printf("  %li steps\n  %li rejected steps\n\n", sys.get_nstep(), sys.get_nrej());
    sys.print();

    //----------------------
    //attempt classification

    printf("\nThe system appears to be...");
    int c = sys.has_oscillation();
    switch (c) {
        case -1:
            printf("stable with an eroding sill.\n\n");
            break;
        case 0:
            printf("oscillating.\n\n");
            break;
        case 1:
            printf("stable after cutoff and  desiccation.\n\n");
            break;
    }

    //---------------------------
    //attempt index approximation

    int rsuc;
    double rzs, rzm, index;
    rsuc = sys.has_root(-250.0, -0.5, -1000.0, -0.5, &rzs, &rzm);
    if ( rsuc == 0 ) {
        printf("fixed point found at\n  zs = %g\n  zm = %g\n", rzs, rzm);
        index = sys.index(rzs, rzm);
        printf("index of fixed point: %g\n\n", index);
    } else {
        printf("no fixed point found\n\n");
    }

    return(0);
}
