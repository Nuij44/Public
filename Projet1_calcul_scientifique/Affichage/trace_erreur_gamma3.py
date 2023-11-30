import numpy as np
import matplotlib.pyplot as plt

res=np.array([55.705087, 28.129810, 23.189416, 22.537224, 22.826640, 23.296657, 23.760265, 24.173838, 24.532532, 24.841784, 25.108997, 25.341117, 25.544046, 25.722628, 25.880793, 26.021725, 26.148012, 26.261768, 26.364731, 26.458342])

plt.plot(np.arange(5, 105 , 5) ,res ,label='erreur L2 sur Gamma3 en fonction de alpha')
plt.title("Tracé du graphe de l'erreur L2 sur Gamma3 en fonction du coefficent de Lavrentier")
plt.legend()
plt.show()
