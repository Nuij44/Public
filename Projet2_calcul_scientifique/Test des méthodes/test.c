#include "outil.h"
#include "fredholm.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

void main(){

    printf("On approche la solution du probleme de l'equation integrale de Fredholm de seconde espece.\n");
    int it_max = 3;
    printf("On approche les series par  les %d premiers termes.\n", it_max);

    //ce bloc implémente la méthode de décomposition de Adomain
    int n = 1000;
    printf("On calcul jusqu'au %d ieme element de la suite.\n", n);
    double alpha = 1;
    double H = 1;
    double L = 1;
    double k_test(int it_max, double H, double L, double x, double t){
        return t;
    }
    double f_test(double x){
        return ( exp(x) - x);
    }


    //test avec exemple analytique
    //printf("U_app_adomain(1) = %lf\n", adomain(1, 1, 20, 0, 1, k_test, f_test, H, L, it_max));
    //printf("U_exacte(1) = %lf\n", exp(1));
    //printf("U_app_iteres(1) = %lf\n", noyaux_iteres(1, 1, 2, 0, 1, k_test, f_test, H, L, it_max));


    FILE *f;
    FILE *g;
    f = fopen("test_adomain.txt", "w");
    g = fopen("test_noyaux_iteres.txt", "w");

    double result;
    int i;
    //on vérifie la bonne ouverture du fichier.
    if (f==NULL){
        printf("Erreur lors de l'ouverture du fichier test_adomain.txt");
        exit(1);
    }
    if (g==NULL)
    {
        printf("Erreur lors de l'ouverture du fichier test_noyaux iteres.txt");
        exit(1);
    }
    

    //on écrit dans les fichiers le resultat approché selon la méthode demandée.
    double res;
    for ( i = 0; i <= 20; i++)
    {
        res = adomain(-1 + i*0.1, -1 + i*0.1, n, 0, 1, k_test, f_test, H, L, it_max);
        fprintf(f,"%lf, ", res);
        res = noyaux_iteres(-1 + i*0.1, -1 + i*0.1, n, 0, 1, k_test, f_test, H, L, it_max);
        fprintf(g,"%lf, ", res);
    }
    fclose(f);
    fclose(g);
}