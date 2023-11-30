#include <stdio.h>
#include <math.h>

//On déclare ici les fonctions associées à la solution excate telles que les conditions aux bords ou la solution exacte

double u_ex(double x, double y){
    double pi = 3.14159265358979323851;
    return (cosh(pi*y)*cos(pi*x));
}

double f0(double x){
    double pi = 3.14159265358979323851;
    return (cos(pi*x));
}

double q0(double x){
    return 0;
}

double f3(double x, double lavrentier){ // On ajoute un paramètre en entrée pour pouvoir injecter la solution exacte pour f_3 et ainsi vérifier l'algorithme de calcul de T
    double H = 2.0;
    double pi = 3.14159265358979323851;
    return (cosh(pi*H)*cos(pi*x) + 0*lavrentier);
}
