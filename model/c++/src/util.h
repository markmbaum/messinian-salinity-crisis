#ifndef UTIL_H_
#define UTIL_H_

//! \file util.h

#include <cmath>
#include <string>
#include <vector>
#include <cstdio>
#include <iostream>

//!max of two numbers
double max (double a, double b);
//!min of two numbers
double min (double a, double b);

//!checks if two numbers are very close to each other
bool is_close (double a, double b, double thresh=1e-4);

//!creates an evenly spaced vector of values over a range
std::vector<double> linspace (double a, double b, long n);

//!create an logarithmically spaced vector of values over a range
std::vector<double> logspace (double a, double b, long n);

//!checks if a file can be written to
void check_file_write (const char *fn);

//!writes a float64 vector to a binary file
void write_double (const char *fn, double *a, long size);

//!writes a float64 array to a binary file
void write_double (const std::string &fn, std::vector<double> a);

//!writes an int array to a binary file
void write_int (const char *fn, int *a, long size);

#endif
