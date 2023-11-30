import numpy as np
import matplotlib.pyplot as plt

res=np.array([66.631267, 32.007921, 24.347715, 22.605743, 22.473317, 22.783322, 23.193481, 23.599140, 23.969834, 24.299571, 24.590413, 24.846752, 25.073271, 25.274259, 25.453426, 25.613912, 25.758341, 25.888908, 26.007450, 26.115507])

plt.plot(np.arange(5, 105 , 5) ,res ,label='erreur L2 sur Gamma3 en fonction de alpha')
plt.title("Tracé du graphe de l'erreur L2 sur Gamma3 en fonction du coefficent de Lavrentier")
plt.legend()
plt.show()
