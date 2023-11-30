*******************************************

Le dossier "Affichage" contient trois codes pyhton qui permettent de tracer les figures utilisées dans le rapport.

-"trace_erreur_gamma3.py" permet de tracer ||u_ex - T^alpha|| sur Gamma_3 en fonction de alpha. Pour cela il faut copier le contenu du fichier "erreur_sur_Gamma3.txt" dans le vecteur nommé 'res'.

-"trace_solution.py" trace sur Omega la solution exacte, la solution approchée ou la différence entre les deux au choix de l'utilisateur. Ce choix ce fait en mettant à un la variable 'solex' pour la solution exacte, 'solapp' pour la solution approchée ou 'diffsol' pour la différence entre les solutions.
Pour le bon fonctionnement le contenu de "solution_exacte_omega.txt" doit être mis dans la variable 'exacte', "solution_calc_omega.txt" dans 'calc' et "difference_calc_exacte.txt" dans 'diff'.
La figure obtenue est un nuage de point n'ayant pas réussi à construire une surface avec le module "plot_suface()".

-"trace_frontiere_gamma.py" trace la comparaison entre la frontière exacte et la frontière approchée. Pour cela dans le vecteur 'Gamma' il faut copier le contenu de "frontiere_Gamma.txt". De plus, il faut définir la fonction 'Gamma_ex(x)' avec l'expression analytique de phi(x) comme décrit dans le rapport.

*******************************************

Le dossier "rapport" contient juste le code latex, les images et le pdf du rapport.

*******************************************

Enfin le dossier "Code_C" contient les fichiers en .c, .txt et .h.

-"outil.c" contient la quadrature de Gauss-Legendre faites sur 5 points avec plusieurs formes suivant que l'intégrande possède en plus un paramètre entier tel que dans des séries ou un paramètre réel tel que le coefficeint de Lavrentier.
Il y a aussi l'algorithme de Newton pour résoudre le problème non linéaire de la détermination de la frontière Gamma.

-"solex.c" sont les données initiales q0, f0 et la solution exacte u_ex ainsi que f3 exacte servant à valider certains blocs de code lors de sa constructions.

-"coefficients.c" définit les coefficients de Fourrier dans le calcul de la solution approchée. On a définit séparément A0 et B0 des itérés suivant pour pouvoir initialiser efficacement le calcul de la solution approchée.

-"frontière" permet de choisir la forme exacte de la frontière Gamma en explicitant la fonction phi.

-"main" contient l'algorithme principal qui permet de choisir la finesse du maillage de Omega avec N_L, le alpha choisit avec alpha_op et la limite du nombre de terme calculés dans les séries avec itmax_serie.

-Enfin les différents fichiers textes contiennent des sauvegardes des résultats calculés pour pouvoir ensuite tracer les figures associées dans le rapport.