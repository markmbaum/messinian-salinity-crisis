//! \file main_sweep.cc

#include <string>
#include <vector>
#include <cstdio>
#include <cstdlib>

#include "omp.h"

#include "util.h"
#include "msc_gcv.h"

//!runs a sweep over ranges of key parameters, searching for oscillatory results
/*!
+ sweeps over `kb`, `tauc`, `Cw`, `U`, `a`, and `L`
+ requires two command line input arguments, the number of values for each of the varied parameters and the output directory
+ each parameter set is tested with the `has_oscillation` method of a `MscGcv` object
*/
int main (int argc, char **argv) {

    //output parameter table
    std::string dirout, fnout;
    FILE *ofile;
    //indices and counters
    long unsigned i, j, temp, idx, nrow;
    //classification vector
    int *cla;
    //parameter table
    double **ptab;
    //get number of threads being used
    int nthread = omp_get_max_threads();
    //array of integrators
    MscGcv *sys = new MscGcv[nthread];

    //--------------------------------------------------------------------------
    // IMPORTANT INPUT VARIABLES

    //output directory, if provided
    if (argc != 3) {
        printf("\ninvalid number of command line arguments\n");
        printf("must supply:\n  1) number of values for each param\n  2) output directory\n");
        exit(EXIT_FAILURE);
    }
    int N = std::stoi(argv[1]);
    dirout = argv[2];
    printf("output directory = '%s'\n", dirout.c_str());

    //number of parameters to sweep over (must be length of pname and pvec)
    int nparam = 6;
    //parameter names
    std::vector< std::string > pname = {"kb", "tauc", "Cw", "U", "a", "L"};
    //parameter values
    std::vector< std::vector<double> > pvec = {
        logspace(-8, -4, N), // kb [m/yr*Pa^-a]
        linspace(25, 100, N), // tauc [Pa]
        linspace(0.5, 10, N), // Cw [?]
        linspace(0.49, 4.9, N), // U [mm/yr
        linspace(1.0, 2.0, N), // a [-]
        linspace(50e3, 150e3, N), // L [m]
    };
    //values of parameters that aren't being varied
    for (int i=0; i<nthread; i++) {
        //set non-varying parameters
        sys[i].n = 0.05;
        sys[i].P = 0.6/YRSEC;
        sys[i].E = 1.2/YRSEC;
        sys[i].R = 4500.0 + 12000.0;
    }

    //print a parameter range summary
    printf("\nparameter range summary:\n\n");
    printf("            | min    | max    | count\n");
    printf("------------|--------|--------|-------\n");
    for (int i=0; i<nparam; i++) {
        printf(" %10s | %6g | %6g | %5d\n",
            pname[i].c_str(),
            pvec[i][0],
            pvec[i][pvec[i].size()-1],
            (int)pvec[i].size()
        );
    }

    //write parameter values
    fnout = dirout + "/parameters.txt";
    check_file_write(fnout.c_str());
    ofile = fopen(fnout.c_str(), "w");
    for (i=0; i<pname.size(); i++) {
        if (i != 0)
            fprintf(ofile, "\n");
        fprintf(ofile, "%-16s", pname[i].c_str());
        for (j=0; j<pvec[i].size(); j++) {
            fprintf(ofile, " %-16g", pvec[i][j]);
        }
    }
    fprintf(ofile, "\n");
    fclose(ofile);
    printf("\nparameter values written to: %s\n", fnout.c_str());

    //--------------------------------------------------------------------------
    // PARAMETER TABLE

    //compute total number of rows
    nrow = 1;
    for (i=0; i<pvec.size(); i++)
        nrow *= pvec[i].size();
    //allocate table columns
    ptab = new double*[nrow];
    for (i=0; i<nrow; i++)
        ptab[i] = new double[nparam];
    //fill in all combinations
    for (i=0; i<nrow; i++) {
        temp = i;
        for (j=0; j<pvec.size(); j++) {
            //find proper index in the parameter vectors
            idx = temp % pvec[j].size();
            temp /= pvec[j].size();
            //set parameter value
            ptab[i][j] = pvec[j][idx];
        }
    }

    //--------------------------------------------------------------------------
    // INTEGRATIONS

    //results array
    cla = new int[nrow];
    //classify in parallel
    printf("\nstarting %lu classifications with %d threads\n", nrow, nthread);
    double ts = omp_get_wtime();
    #pragma omp parallel for schedule(dynamic)
    for (i=0; i<nrow; i++) {
        //get the thread number
        int n = omp_get_thread_num();
        //set varying parameters
        sys[n].kb =         ptab[i][0]/YRSEC;
        sys[n].tauc =       ptab[i][1];
        sys[n].Cw =         ptab[i][2];
        sys[n].U =          ptab[i][3]/MMYR;
        sys[n].a =          ptab[i][4];
        sys[n].L =          ptab[i][5];
        //classify
        cla[i] = sys[n].has_oscillation();
    }
    delete [] sys;
    ts = omp_get_wtime() - ts;
    printf("classifications finished\n  %g seconds total\n  %g seconds per trial",
        ts, ts/nrow);

    //--------------------------------------------------------------------------
    // OUTPUT

    //print basic classification stats
    long unsigned nsta = 0,
                  nosc = 0,
                  ndes = 0;
    for (i=0; i<nrow; i++) {
        switch (cla[i]) {
          case -1: nsta++; break;
          case  0: nosc++; break;
          case  1: ndes++; break;
        }
    }

    printf("\nclassification summary:\n");
    printf("  %g %% stable with eroding sill\n", 100*double(nsta)/nrow);
    printf("  %g %% oscillating\n", 100*double(nosc)/nrow);
    printf("  %g %% stable after cutoff and desiccation\n", 100*double(ndes)/nrow);

    //write output table and parameter vectors
    fnout = dirout + "/trials.csv";
    check_file_write(fnout.c_str());
    ofile = fopen(fnout.c_str(), "w");
    fprintf(ofile, "trial");
    for (j=0; j<pname.size(); j++) fprintf(ofile, ",%s", pname[j].c_str());
    fprintf(ofile, ",classification\n");
    for (i=0; i<nrow; i++) {
        fprintf(ofile, "%lu", i);
        for (j=0; j<pvec.size(); j++) fprintf(ofile, ",%g", ptab[i][j]);
        fprintf(ofile, ",%d", cla[i]);
        fprintf(ofile, "\n");
    }
    //close the file
    fclose(ofile);
    printf("\ntable written to: %s\n", fnout.c_str());
    printf("   1 = stable solution after cutoff and dessication\n");
    printf("   0 = oscillating solution\n");
    printf("  -1 = stable solution with an eroding sill\n");

    //--------------------------------------------------------------------------

    for (i=0; i<nrow; i++) delete [] ptab[i];
    delete [] ptab;
    delete [] cla;

    return(0);
}
