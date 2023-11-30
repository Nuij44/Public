#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//fonction d'approximation des intégrales
double integrate(double borne_inf, double borne_sup, double (f)(double)){
    // on utilise la méthode de quadrature de Gauss-Legendre avec 5 points

    double xi[] = {-(1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), 0, (1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0))};
    double wi[] = {(322.0 - 13.0*sqrt(70.0))/900.0, (322.0 + 13.0*sqrt(70.0))/900.0, 128.0/225.0, (322.0 + 13.0*sqrt(70.0))/900.0, (322.0 - 13.0*sqrt(70.0))/900.0};
    int taille_i = 5;
    double result = 0;
    int i=0;
    for ( i = 0; i < taille_i; i++)
    {
        result = result + wi[i]*f(((borne_sup - borne_inf)/2)*xi[i] + ((borne_sup + borne_inf)/2));
    }
    result = ((borne_sup - borne_inf)/2)*result;
    return result;

}

double integrate_m(double borne_inf, double borne_sup, int m, double (f)(double, int)){
    // on utilise la méthode de quadrature de Gauss-Legendre avec 5 points

    double xi[] = {-(1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), 0, (1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0))};
    double wi[] = {(322.0 - 13.0*sqrt(70.0))/900.0, (322.0 + 13.0*sqrt(70.0))/900.0, 128.0/225.0, (322.0 + 13.0*sqrt(70.0))/900.0, (322.0 - 13.0*sqrt(70.0))/900.0};
    int taille_i = 5;
    double result = 0;
    int i=0;
    for ( i = 0; i < 5; i++)
    {
        result = result + wi[i]*f(((borne_sup - borne_inf)/2)*xi[i] + ((borne_sup + borne_inf)/2), m);
    }
    result = ((borne_sup - borne_inf)/2)*result;
    return result;
}

double integrate_alpha(double borne_inf, double borne_sup, double lavrentier, double (f)(double, double)){
    // on utilise la méthode de quadrature de Gauss-Legendre avec 5 points

    double xi[] = {-(1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), 0, (1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0))};
    double wi[] = {(322.0 - 13.0*sqrt(70.0))/900.0, (322.0 + 13.0*sqrt(70.0))/900.0, 128.0/225.0, (322.0 + 13.0*sqrt(70.0))/900.0, (322.0 - 13.0*sqrt(70.0))/900.0};
    int taille_i = 5;
    double result = 0;
    int i=0;
    for ( i = 0; i < 5; i++)
    {
        result = result + wi[i]*f(((borne_sup - borne_inf)/2)*xi[i] + ((borne_sup + borne_inf)/2), lavrentier);
    }
    result = ((borne_sup - borne_inf)/2)*result;
    return result;


}
//On implémente l'algorithme de newton
 double Newton(int max, double H, int m, double (f)(double, int), double (f_prime)(double, int)){
    double U_k;
    int i;
    //on initialise la suite avec un nombre aléatoire compris entre 0 et L
    /*
    srand(time(NULL));
    U_k = H*rand()/(RAND_MAX+1.0);
    */
   U_k = 0.5;
    //On itère jusqu'à max pour éviter d'avoir à faire de test logique
    for ( i = 0; i < max; i++)
    {
        U_k = U_k - (f(U_k, m)/f_prime(U_k, m));
    }
    /*
    if ((U_k < 0.0) || (U_k > H))
    {
        U_k = Newton(max, H, m, f, f_prime);
    }
    */
    
    return U_k;
 }