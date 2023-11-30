#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "outil.h"
#include "coefficients.h"
#include "solex.h"
#include "frontiere.h"
#include "fredholm.h"
#include <stdbool.h>

void main(){
    int choix;
    double L = 1.0;
    double H = 1.0;
    double pi = 3.14159265358979323851;
    double alpha_op = 21 ;
    int itmax_serie = 20; 
    int N_L = 10;
    int log = 0; // mettre 1 pour afficher les logs et 0 sinon
    int log_series = 0; // mettre 1 pour afficher les logs dans le calcul des séries et 0 sinon
    int ecrire_erreur_gamma3 = 0 ; //Mettre 1 pour afficher l'erreur sur le segment Gamma 3
    int recherch_alpha = 0; //Mettre 1 pout remplir le fichier pour different alpha l'erreur commise sur Gamma 3
    int solution_Omega = 0; //Mettre 1 pour remplir les fichiers des valeurs sur les points du maillage de Omega
    int frontiere_gamma = 1; //Mettre 1 pour passer à la recherche de la frontiere Gamma et penser à rentrer T0

    printf("On calcul la solution approchee pour le probleme geometrique inverse avec H = %lf et L = %lf\n",H,L);
    printf("Veuillez choisir ce que l'algorithme recherche en entrant votre choix dans le code source:\n");
    printf("1- Afficher l'erreur commise par rapport a la solution exacte sur le Gamma3.\n");
    printf("2- Remplir le fichier 'erreur_sur_Gamma3.txt' de la norme de l'erreur en fonction de differents coefficients de Lavrentier.\n");
    printf("3- Calculer la solution exacte, la solution approchee et la difference des deux sur le maillage et sauvegarder cela dans trois fichiers.\n");
    printf("4- Recherche de la frontiere Gamma.\n");
    
    printf("Choix de la méthode d'approximation pour f3(x) à changer dans le code source :\n");
    printf("1- Méthode de décomposition d'Adomain\n");
    printf("2- Méthode des noyaux itérés\n");
    
    int Methode_Adomain = 0;    //Mettre 1 pour approcher f3 par une décomposition d'Adomain
    int Methode_noyaux_iteres = 1;  //Mettre 1 pour approcher f3 par la méthode des noyaux itérés
    int nb_iteration_methode = 10;  //Nombre d'iteration de la méthode chosit pour la résolution de l'équation de Fredholm

    printf("On résoud l'équation de Fredholm avec %d itérations\n", nb_iteration_methode);

    //On définit h
    double h(double xh){
        double somme = 0.0;
        int k;
        //on constuit l'intégrande
        double integrande_h(double t, int i){
            return (cos(i*pi*t/L)*f0(t));
        }
        for (k = 1; k < itmax_serie; k++)
        {
            //printf("terme de la somme = %lf et  k = %d\n", ((k*cos(xh*k*pi/L))/sinh(k*pi*H/L))*integrate_m(0, L, k, integrande_h), k);
            somme = somme + ((k*cos(xh*k*pi/L))/sinh(k*pi*H/L))*cosh(H*k*pi/L)*integrate_m(0, L, k, integrande_h);
            if (log_series == 1)
        {
            printf("Série pour le calcul de h = %lf\n", somme);
        }
        }
        return (-q0(xh) + (1/H)*A_0(L) + (2*pi/(L*L))*somme);
    }
    
    //On définit k(x,t)
    double k(int it_max, double H, double L, double x, double t){
        int i;
        double somme = 0;
        double pi = 3.14159265358979323851;
        for ( i = 1; i < it_max; i++)
        {
            somme = somme + (i*cos(i*pi*x/L)*cos(i*pi*t/L))/sinh(i*pi*H/L);
        }
        return (1/(H*L) + 2*pi/(L*L)*somme);
    }

    //On construit une approximation de f3(x) 
    double f3_app(double x, double lavrentier){

        //On consruit maintenant les éléments notés f(x) et lambda dans le rapport de Mr. Bounefla Chemes Eddine
        double f_rapport(double xf){
            return (h(xf)/lavrentier);
        }

        double lambda = -1/lavrentier;

        //On fait un arbre de décision selon la méthode choisit
        if (Methode_Adomain == 1)
        {
            return adomain( x, lambda, nb_iteration_methode, 0, L, k, f_rapport, H, L, itmax_serie);
        }
        else if (Methode_noyaux_iteres == 1)
        {
            return noyaux_iteres(x, lambda, nb_iteration_methode, 0, L, k, f_rapport, H, L, itmax_serie);
        }
    }



    //on construit T^alpha(xi) pour pouvoir écrire (xi, t^alpha(xi))
    double T_alpha(double xi, double yi, double lavrentier, int itmax){
        double T = A_0(L) + B_0(H, L, lavrentier, f3_app)*yi;
        int m;
        for ( m = 1; m < itmax; m++)
        {
            T = T + (A_m(H, L, m, lavrentier, f3_app)*exp(yi*m*pi/L) + B_m(H, L, m, lavrentier, f3_app)*exp(-yi*m*pi/L))*cos(xi*pi*m/L);
            if (log_series == 1)
            {
                printf(" A%d = %lf", m, A_m(H, L, 1, lavrentier, f3_app));
                printf(" B%d = %lf", m, B_m(H, L, 1, lavrentier, f3_app));
                printf("T_%d = %lf", m, T);
            }
            
            
        }
        return T;
    }
    
    //On écrit les (xi, t^alpha(xi)) dans le fichier resultats.txt
    
    FILE *f;
    f = fopen("resultats.txt", "w");
    
    double result;
    int i;
    //on vérifie la bonne ouverture du fichier.
    if (f==NULL){
        printf("Erreur lors de l'ouverture d'un fichier");
        exit(1);
    }

    //on écrit dans resultats.txt le couple demandé.
    for ( i = 0; i <= N_L; i++)
    {
        result = T_alpha(i*(L/N_L), H, alpha_op, itmax_serie);
        fprintf(f, "(%f, %lf)       (%f, %lf)\n", (i*(L/N_L)), result, (i*(L/N_L)), u_ex((i*(L/N_L)),H));
    }
    fclose(f);

    //On calcul la norme L2 de l'erreur sur Gamma3
    double norme_L2_G3(double lavrentier){
        double norme = 0;
        double erreur;
        int i;
        for ( i = 0; i <= N_L; i++)
        {
            erreur = u_ex(i*(L/N_L), H) - T_alpha(i*(L/N_L), H, lavrentier, itmax_serie);
            norme = norme + erreur*erreur;
        }
        return  sqrt(norme);
    }
    if (ecrire_erreur_gamma3 == 1)
    {
        printf("||u_ex - T_alpha||L2(G3) = %lf\n", norme_L2_G3(alpha_op));
        if (log==1)
        {
            printf("T_exact (0.3, H) = %lf\n", u_ex(0.3, 0));
            printf("T approx = %lf\n", T_alpha(0.3, 0, alpha_op, itmax_serie));
            printf("première approx m=0, t_calc = %lf\n", ((2/pi) + 2.0*(cosh(2*pi) - 1)/pi));
        }
    }
    
    


    
    if (recherch_alpha = 1)
    {
        FILE *g;
        g = fopen("erreur_sur_Gamma3.txt", "w");
    

        //on vérifie la bonne ouverture du fichier.
        if (g==NULL){
            printf("Erreur lors de l'ouverture d'un fichier");
            exit(1);
        }
        //double coefficient_lavrentier[] ={pow(10,-1), 2*pow(10,-1), 3*pow(10,-1), 4*pow(10,-1), 5*pow(10,-1), 6*pow(10,-1), 7*pow(10,-1), 8*pow(10,-1), 9*pow(10,-1), 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
        double coefficient_lavrentier[] = {5, 10 ,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
        double norme_er = 0;
        for ( i = 0; i <= 19; i++)
        {
            norme_er = norme_L2_G3(coefficient_lavrentier[i]);
            fprintf(g,"%lf, ", norme_er);
        }
        fclose(g);
    }
    


    //On calcul maintenant la solution sur un maillage de Omega N_LxN_L pour la solution exacte, la solution approchée et leur différence absolue
    if (solution_Omega = 1)
    {
        FILE *calc;
        FILE *diff;
        diff = fopen("difference_calc_exacte.txt", "w");
        calc = fopen("solution_calc_omega.txt", "w");
    

        //on vérifie la bonne ouverture du fichier.
        if (calc==NULL){
            printf("Erreur lors de l'ouverture d'un fichier");
            exit(1);
        }
        if (diff==NULL)
        {
            printf("Erreur lors de l'ouverture d'un fichier");
            exit(1);
        }
        
        
        int j;
        double stockage_calc;
        double stock;
        for ( i = 1; i <= N_L; i++)
        {
            for ( j = 1; j <= N_L; j++)
            {
                stockage_calc = T_alpha(L*j/N_L, H*i/N_L, alpha_op, itmax_serie);
                fprintf(calc,"%lf, ", stockage_calc);
                stock = fabs(stockage_calc - u_ex(L*j/N_L, H*i/N_L));
                fprintf(diff, "%lf, ", stock);
            }
        }
        fclose(calc);
        fclose(diff);

        FILE *exacte;
        exacte = fopen("solution_exacte_omega.txt", "w");
    

        //on vérifie la bonne ouverture du fichier.
        if (exacte==NULL){
            printf("Erreur lors de l'ouverture d'un fichier");
            exit(1);
        }
        for ( i = 1; i <= N_L; i++)
        {
            for ( j = 1; j <= N_L; j++)
            {
                fprintf(exacte,"%lf, ", u_ex(L*j/N_L, H*i/N_L));
            }
        }
        fclose(exacte);
    }
    

    //On passe ici sur la recherche de la frontiere Gamma
    if (frontiere_gamma == 1)
    {
        double f_newton(double x, double y){
            return (T_alpha(x, y, alpha_op, itmax_serie));
        }
        //On a calculé analytiquement la dérivée 
        double f_newton_prime(double x, double y){
            double s = B_0(H, L, alpha_op, f3_app);
            for ( i = 1; i <= itmax_serie; i++)
            {
                s = s + cos(x*i*pi/L)*(i*pi/L)*(A_m(H, L, i, alpha_op, f3_app)*exp(i*pi*y/L) - B_m(H, L, i, alpha_op, f3_app)*exp(-i*pi*y/L));
            }
            return s;
        }

        //on boucle sur la discretisation de [0,L] avec l'algorithme de Newton
        double Gamma[N_L + 1];

        double F(double y, int k){
            return (T_alpha(k*L/N_L, y, alpha_op, itmax_serie) - T0(k*L/N_L));
        }

        double F_prime(double y, int k){
            return (f_newton_prime(k*L/N_L, y));
        }
        int ii;
        for ( ii = 0; ii <= N_L; ii++)
        {
            Gamma[ii] = Newton(10, L, ii, F, F_prime);
            printf("Gamma en  x = %lf et y = %lf", ii*L/N_L, Gamma[ii]);
        }
        
        //On écrit le résultat dans un fichier
        FILE *front;
        front = fopen("frontiere_Gamma.txt", "w");
    

        //on vérifie la bonne ouverture du fichier.
        if (front==NULL){
            printf("Erreur lors de l'ouverture d'un fichier");
            exit(1);
        }
        for ( i = 0; i <= N_L; i++)
        {
            fprintf(front,"%lf, ", Gamma[i]);
        }
        fclose(front);
    }
    
}