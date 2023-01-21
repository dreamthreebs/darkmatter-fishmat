from matplotlib import rc # latex in matplotlib
from matplotlib import pyplot as plt
import numpy as np
rc('text',usetex=True)
x=np.arange(10)
y=np.geomspace(1,10,10)
plt.plot(x,y)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\chi\chi\rightarrow\gamma\gamma$')
plt.title(r'fsky=0.5,$\left\langle \sigma v\right\rangle/ m_\chi \left[cm^3s^{-1}GeV^{-1}\right]$')
plt.show()

