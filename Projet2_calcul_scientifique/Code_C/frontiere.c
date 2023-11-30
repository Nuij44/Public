#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "outil.h"
#include "coefficients.h"
#include "solex.h"
#include "frontiere.h"
#include <stdbool.h>


//T0 est à choisir pour tester l'algorithme
double T0(double x){
    double pi = 3.14159265358979323851;
    double phi(double x){
        return (0.2*sin(pi*x) + 0.2);
    }
    return (u_ex(x, phi(x)));
}




