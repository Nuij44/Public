#ifndef fredholm_h_
#define fredholm_h_

double adomain(double x, double lambda, int n, double a, double b, double (k)(int, double, double, double, double), double (f)(double), double H, double L, double it_max);
double noyaux_iteres(double x, double lambda, int n, double a, double b, double (k)(int, double, double, double, double), double (f)(double), double H, double L, double it_max);
    

#endif