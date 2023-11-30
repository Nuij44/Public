import numpy as np
import matplotlib.pyplot as plt

adomain=np.array([0.590101, 0.622777, 0.655043, 0.687573, 0.721119, 0.756531, 0.794764, 0.836905, 0.884185, 0.938012, 1.000000, 1.072013, 1.156218, 1.255153, 1.371825, 1.509832, 1.673547, 1.868368, 2.101096, 2.380512, 2.718281])
noyaux=np.array([1.9 ,1.812325, 1.766340, 1.725519, 1.690412, 1.661632, 1.639864, 1.625876, 1.620528, 1.624791, 1.639758, 1.666667, 1.706925, 1.762143, 1.834173, 1.925158, 2.037610, 2.174500, 2.339394, 2.536652, 2.771724 ])
X=np.arange(-1,1.1,0.1)
plt.plot(X ,adomain ,label='Méthode de decomposition de Adomain')
plt.plot(X, noyaux, label='Méthode des noyaux itérés')
plt.plot(X, np.exp(X), label='Solution exacte')
plt.title("Comparaison des méthodes sur le cas test.")
plt.legend()
plt.show()
