import numpy as np
import matplotlib.pyplot as plt

pi = 3.14159265358979323851
Gamma = np.array([0.552817, 0.598982, 0.633340, 0.658888, 0.681923, -0.192421, 0.681923, 0.658888, 0.633340, 0.598982, 0.552817, ])
def Gamma_ex(x):
    return (0.2 + 0.5*np.sin(np.pi*x))
plt.ylim(0,1)
plt.xlim(0,1)
plt.xlabel('x')
plt.ylabel('y')
parameter=1/(np.size(Gamma) - 1)
axe_x = np.arange(0, 1 + parameter, parameter)
plt.plot(axe_x ,Gamma ,label='Frontière approchée Gamma')
plt.plot(axe_x ,Gamma_ex(axe_x) ,label='Frontière exacte Gamma')
plt.title("Approximation de la frontière Gamma dans Omega par les noyaux itérés")
plt.legend()
plt.show()
