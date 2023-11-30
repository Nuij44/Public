//header des coefficients

#ifndef coefficients_h_
#define coefficients_h_

double A_0(double L);
double B_0(double H, double L, double lavrentier, double (f_3)(double, double));
double A_m(double H, double L, int m, double lavrentier, double (f_3)(double, double));
double B_m(double H, double L, int m, double lavrentier, double (f_3)(double, double));

#endif