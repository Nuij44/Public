#include <stdio.h>
#include <math.h>
#include "fredholm.h"
#include "outil.h"



double adomain(double x, double lambda, int n, double a, double b, double (k)(int, double, double, double, double), double (f)(double), double H, double L, double it_max){
    int i;
    int j;
    int l;
    //On définit les points xi de Legendre
    double ti[] = {-(1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), 0, (1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0))};
    
    //On renormalise pour être sur l'intervalle [0,1]
    for ( i = 0; i < 5; i++)
    {
        ti[i] = ((b - a)/2)*ti[i] + (a + b)/2;
    }

    //On initialise la somme de U_n
    double somme = f(x);
    //printf("somme partielle = %lf\n", somme);

    //on test le cas où f(x)=0
    if (somme == 0)
    {
        somme = 1e-5;
    }

    //On initialise la boucle sur n
    double U[5];
    double U_plus[5];
    double W[5];
    double Z[5];
    for ( i = 0; i < 5; i++)
    {
        for ( l = 0; l < 5; l++)
            {
                W[l] = k(it_max, H, L, ti[i], ti[l])*f(ti[l]);
            }
        U[i] = lambda*integrate_vect(a,b,W);
        //printf("U[%d] = %lf\n", i, U[i]);
        //printf("k(x,ti[%d]) = %lf\n", i, k(it_max, H, L, x, ti[i]));
    }

    for ( i = 0; i < 5; i++)
    {
        Z[i] = k(it_max, H, L, x, ti[i])*f(ti[i]);
    }
    somme = somme + lambda*integrate_vect(a,b,Z);

    //On fait la boucle sur n
    for ( i = 1; i < n; i++)
    {
        for ( j = 0; j < 5; j++)
        {
            for ( l = 0; l < 5; l++)
            {
                W[l] = k(it_max, H, L, ti[j],ti[l])*U[l];
            }
            
            U_plus[j] = lambda*integrate_vect( a, b, W);
            //printf("U_plus[%d] = %lf\n", j, U_plus[j]);
        }

        for ( l = 0; l < 5; l++)
        {
            W[l] = k(it_max, H, L, x, ti[l])*U_plus[l];
        }
        
        somme = somme + lambda*integrate_vect( a, b, W);
        //printf("somme partielle = %lf\n", somme);
        //printf(" n = %d\n", i);
        
        //On fait une seconde boucle pour décaler les indices
        for ( j = 0; j < 5; j++)
        {
            U[j] = U_plus[j];
        }
    }
    return somme;
}

double noyaux_iteres(double x, double lambda, int n, double a, double b, double (k)(int, double, double, double, double), double (f)(double), double H, double L, double it_max){
    
    //On utilise quasiment le même code que pour adomain mais on évalue en x seulement le dernier itéré Un

    int i;
    int j;
    int l;
    //On définit les points xi de Legendre
    double ti[] = {-(1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0)), -(1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), 0, (1.0/3.0)*sqrt(5 - 2*sqrt(10.0/7.0)), (1.0/3.0)*sqrt(5 + 2*sqrt(10.0/7.0))};
    
    //On renormalise pour être sur l'intervalle [a,b]
    for ( i = 0; i < 5; i++)
    {
        ti[i] = ((b - a)/2)*ti[i] + (a + b)/2;
    }

    //On initialise la boucle sur n
    double U[5];
    double U_plus[5];
    double W[5];
    double resultat;
    for ( i = 0; i < 5; i++)
    {
        for ( l = 0; l < 5; l++)
            {
                W[l] = k(it_max, H, L, ti[i], ti[l])*f(ti[l]);
            }
        U[i] = f(ti[i]) + lambda*integrate_vect(a,b,W);
        //printf("U[%d] = %lf\n", i, U[i]);
        //printf("k(x,ti[%d]) = %lf\n", i, k(it_max, H, L, x, ti[i]));
        
    }
    for ( j = 0; j < 5; j++)
        {
            W[j] = k(it_max, H, L, x, ti[j])*U[j];
        }
    
    
    resultat = f(x) + integrate_vect(a,b,W);
    //printf("U_app_1(1) = %lf\n", resultat);

    //On fait la boucle sur n
    for ( i = 2; i < (n - 1); i++)
    {
        for ( j = 0; j < 5; j++)
        {
            for ( l = 0; l < 5; l++)
            {
                W[l] = k(it_max, H, L, ti[j],ti[l])*U[l];
            }
            
            U_plus[j] = f(ti[j]) + lambda*integrate_vect( a, b, W);
            //printf("U_plus[%d] = %lf\n", j, U_plus[j]);
        }

        for ( j = 0; j < 5; j++)
        {
            W[j] = k(it_max, H, L, x, ti[j])*U[j];
        }
    
    
        resultat = f(x) + integrate_vect(a,b,W);
       // printf("U_app_%d(1) = %lf\n", i, resultat);

        //On fait une seconde boucle pour décaler les indices
        for ( j = 0; j < 5; j++)
        {
            U[j] = U_plus[j];
        }
    }

    //on évalue l'itéré suivant en x
    for ( i = 0; i < 5; i++)
    {
        W[i] = k(it_max, H, L, x, ti[i])*U[i];
    }
    
    
    resultat = f(x) + integrate_vect(a,b,W);
    return resultat;

}