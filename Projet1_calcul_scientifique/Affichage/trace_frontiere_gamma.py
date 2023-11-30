import numpy as np
import matplotlib.pyplot as plt

pi = 3.14159265358979323851
Gamma = np.array([0.290138, 0.537699, 0.719014, 0.861973, 0.981066, -0.000000, 0.944388, 0.846073, 0.709531, 0.531129, 0.285388, ])
def Gamma_ex(x):
    return (0.2*np.sin(pi*x)+0.5)
plt.ylim(0,1)
plt.xlim(0,1)
plt.xlabel('x')
plt.ylabel('y')
parameter=1/(np.size(Gamma) - 1)
axe_x = np.arange(0, 1 + parameter, parameter)
plt.plot(axe_x ,Gamma ,label='Frontière approchée Gamma')
plt.plot(axe_x ,Gamma_ex(axe_x) ,label='Frontière exacte Gamma')
plt.title("Approximation de la frontière Gamma dans Omega")
plt.legend()
plt.show()
