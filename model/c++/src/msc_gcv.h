//! \file msc_gcv.h
/*!
\mainpage msc_gcv

This code is a simplified, consolidated, C++ implementation of a prior Messinian Salinity Crisis model:
+ [Garcia-Castellanos, D. & Villaseñor, A. Messinian salinity crisis regulated by competing tectonics and erosion at the Gibraltar arc. Nature 480, 359–363 (2011).](http://www.nature.com/articles/nature10651)

The equations and variables of the model are implemented in the `msc_gcv.h` and `msc_gcv.cc` files as a class called MscGcv. The `linalg.cc` and `newton.cc` files contain standard numerical algorithms used by the model. The `util.cc` file has some miscellaneous useful functions. Finally, the `main_single.cc` and `main_sweep.cc` files are two separate driver files. They are compiled by the Makefile into `bin/single.exe` and `bin/sweep.exe`. The `single.exe` program runs a single integration of the model, writes output into the `out` directory, and attempts to find and display the model's fixed point. The results of `single.exe` can be plotted with `scripts/plot_out.py`. The `sweep.exe` file runs, in parallel, a grid of many integrations over ranges of key parameters, testing whether the model oscillates. The model runs on top of ODE integrators from [libode](https://github.com/wordsworthgroup/libode).

Basic steps to compile the code:
1. download and compile [libode](https://github.com/wordsworthgroup/libode)
2. copy/rename the `_config.mk` file `config.mk`, getting rid of the leading underscore
3. edit `config.mk` to have the proper compiler configuration and path to `libode`
4. run `make` at the command line
5. `single.exe` and `sweep.exe` should be present in the `bin` directory

To run `single.exe` and plot the results:
\code{.sh}
  ./bin/single.exe
  cd scripts
  python plot_out.py
  cd ..
\endcode

To run `sweep.exe`:
\code{.sh}
  ./bin/sweep.exe N out
  cd scripts
  python slice_sweep_combos.py
  cd ..
\endcode
replacing N with the number of values to use for each varied parameter.
*/

#ifndef MSC_GCV_H_
#define MSC_GCV_H_

#include <cmath>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>

#include "util.h"
#include "newton.h"

//header file for ODE integrator class
#include "ode_vern_65.h"

//!seconds in an Earth year
#define YRSEC 31557600.0
//!seconds in 1000 Earth years
#define KYRSEC 31557600000.0
//!converts m/s to mm/yr
#define MMYR 31557600000.0
//!pi
#define PI 3.14159265358979323846264338
//!gravitational acceleration [m/s^2]
#define GRAV 9.8
//!water density [kg/m^3]
#define RHO 1000.0

class MscGcv : public OdeVern65, public Newton {

public:

    //!constructs
    MscGcv ();

    //-----------------
    //system parameters

    //!erodability [m/yr*Pa^a]
    double kb;
    //!critical shear stress [Pa]
    double tauc;
    //!constant for Turowski channel width formula [?]
    double Cw;
    //!uplift rate [m/s]
    double U;
    //!erosion exponent [-]
    double a;
    //!stream length [m]
    double L;
    //!roughness coefficient [s/m^(1/3)]
    double n;
    //!precipitation rate [m/s]
    double P;
    //!evaporation rate [m/s]
    double E;
    //!input from non-Gibraltar rivers [m^3/s]
    double R;

    //----------------------------------
    //Mediterranean area and ocean level

    //!computes area of Mediterranean
    double fAm (double zm);
    //!computes ocean level from Mediterranean level
    double fzo (double zm);

    //-------------------------
    //slope for modified system

    //!calculates slope when zm > zs
    double fS (double zs, double zm);

    //--------------
    //system of ODEs

    //!implements the system of ODEs, computing time derivatives
    void ode_fun (double *solin, double *fout);

    //----------------------------------------------
    //finding fixed pts and checking for oscillation

    //!executes Newton's method over a grid of points to find a fixed point
    /*!
    \param[in] zslo lowest sill level in grid
    \param[in] zshi highest sill level in grid
    \param[in] zmlo lowest Med. level in grid
    \param[in] zmhi highest Med. level in grid
    \param[in] nzs number of evenly spaced sill levels in grid
    \param[in] nzm number of evenly spaced Med. levels in grid
    \param[out] rzs sill level of root or NAN
    \param[out] rzm Med. level of root or NAN
    \return success code, zero for success
    */
    int has_root (double zslo, double zshi,
                  double zmlo, double zmhi,
                  double *rzs, double *rzm,
                  int nzs=50, int nzm=50);

    //!integrates for a very long period, checking model state after shorter intervals to determine if it is stable or oscillating
    /*!
    \param[in] zs0 initial sill level
    \param[in] zm0 initial Mediterranean level
    \param[in] tint shorter integration interval duration
    \param[in] tlim total integration time limit
    \param[in] rootrate maximum magnitude of both time derivatives for assuming stationary state
    \return oscillation code: -1=fixed pt with sill below ocean, 0=oscillating, 1=sill above ocean
    */
    int has_oscillation (double zs0=-60.0, double zm0=0.0,
                         double tint=25*KYRSEC, double tlim=1e8*YRSEC,
                         double rootrate=1e-3/MMYR);

    //!computes the index of a point (normally a fixed point) by brute force
    /*!
    + [https://staff.ul.ie/burkem/Teaching/No%20Cycles.pdf](https://staff.ul.ie/burkem/Teaching/No%20Cycles.pdf)
    \param[in] zs sill level of point
    \param[in] zm Mediterranean level of point
    \param[in] r radius around point to integrate
    \return floating point approximation to point's index (should be very nearly an integer)
    */
    double index (double zs, double zm, double r=0.01);

    //!prints zo, zs, zm
    void print ();

private:

    //params for Mediterranean area and ocean level
    const double c1, a1, c2, a2, Ao;

    //implements system of equations for Newton solver by calling ode_fun
    void f_Newton (double *x, double *f) { ode_fun(x, f); }

};

#endif
