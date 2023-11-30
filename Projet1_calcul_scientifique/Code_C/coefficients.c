#include <stdio.h>
#include <math.h>
#include "outil.h"
#include "solex.h"

//On déclare ici les fonctions associer aux coefficients de Fourrier

double A_0(double L){
    return ((1/L)*integrate(0, L, f0));
}

double B_0(double H, double L, double lavrentier, double (f_3)(double, double)){
    double B = integrate_alpha(0, L, lavrentier, f_3) - integrate(0, L, f0);
    return ((1/(H*L))*B);
}

double A_m(double H, double L, int m, double lavrentier, double (f_3)(double, double)){
    double pi = 3.14159265358979323851;
    //on déclare la fonction que l'on va intégrer
    double f_integrandeA(double x){
        return ((f_3(x, lavrentier) - exp(-H*(m*pi/L))*f0(x))*cos(x*m*pi/L));
    }
    return ((1/(L*sinh(H*m*pi/L)))*integrate(0, L, f_integrandeA));
}

double B_m(double H, double L, int m, double lavrentier, double (f_3)(double, double)){
    double pi = 3.14159265358979323851;
    //on déclare la fonction que l'on va intégrer
    double f_integrandeB(double x){
        return ((exp(H*(m*pi/L))*f0(x) - f_3(x, lavrentier))*cos(x*m*pi/L));
    }
    return ((1/(L*sinh(H*m*pi/L)))*integrate(0, L, f_integrandeB));
}
