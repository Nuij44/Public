//header des outils utilisés dans le projet

#ifndef outil_h_
#define outil_h_

//protype de l'approximation de l'intégrale de f entre borne_inf et borne_sup
double integrate(double borne_inf, double borne_sup, double (f)(double));

//protype de l'approximation de l'intégrale de f avec un paramètre m pour les sommes dans h(x) et f_3(x)
double integrate_m(double borne_inf, double borne_sup, int m, double (f)(double, int));

//prototupe pout l'intégrale avec f_3 approchée utilisant ainsi un paramtre un plus : le coefficient de Lavrentier
double integrate_alpha(double borne_inf, double borne_sup, double lavrentier, double (f)(double, double));

//prototype de l'algorithme de Newton
double Newton(int max, double H, int m, double (f)(double, int), double (f_prime)(double, int));
#endif