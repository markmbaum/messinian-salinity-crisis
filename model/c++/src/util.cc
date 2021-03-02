//! \file util.cc

#include "util.h"

double max (double a, double b) {
    if (a > b) return(a);
    return(b);
}

double min (double a, double b) {
    if (a < b) return(a);
    return(b);
}

bool is_close (double a, double b, double thresh) {

    //get magnitudes
    double absa = fabs(a);
    double absb = fabs(b);
    double absd = fabs(a - b);
    //check relative differnence against a threshold
    if ((absd/absa < thresh) && (absd/absb < thresh))
        return(true);
    //otherwise the numbers aren't close
    return(false);
}

std::vector<double> linspace (double a, double b, long n) {

    //avoid division by zero
    if ( n == 1 ) {
        if ( a != b ) {
            printf("FAILURE: cannot linspace with a single value if the range limits are not identical\n");
            exit(EXIT_FAILURE);
        }
        std::vector<double> v(1, a);
        return(v);
    }
    //spacing
    double h = (b - a)/(n - 1);
    //values
    std::vector<double> v;
    for (long i=0; i<n; i++)
        v.push_back( a + i*h );

    return(v);
}

std::vector<double> logspace (double a, double b, long n) {

    //avoid division by zero
    if ( n == 1 ) {
        if ( a != b ) {
            printf("FAILURE: cannot logspace with a single value if the range limits are not identical\n");
            exit(EXIT_FAILURE);
        }
        std::vector<double> v(1, pow(10, a));
        return(v);
    }
    //spacing
    double h = (b - a)/(n - 1);
    //values
    std::vector<double> v;
    for (long i=0; i<n; i++)
        v.push_back( pow(10, a + i*h) );

    return(v);
}

void check_file_write (const char *fn) {
    FILE* ofile;
    ofile = fopen(fn, "w");
    if (ofile == NULL) {
        std::cout << "FAILURE: cannot open file " << fn << std::endl;
        exit(EXIT_FAILURE);
    }
    fclose(ofile);
}

void write_double (const char *fn, double *a, long size) {
    FILE* ofile;
    check_file_write(fn);
    ofile = fopen(fn, "wb");
    fwrite(a, sizeof(double), size, ofile);
    fclose(ofile);
}

void write_double (const std::string &fn, std::vector<double> a) {
    write_double(fn.c_str(), a.data(), long(a.size()));
}

void write_int (const char *fn, int *a, long size) {
    FILE* ofile;
    check_file_write(fn);
    ofile = fopen(fn, "wb");
    fwrite(a, sizeof(int), size, ofile);
    fclose(ofile);
}
