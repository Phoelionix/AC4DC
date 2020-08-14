# Script to generate relevant parameters for Gaussian quadrature integrals

import numpy as np
import numpy.polynomial as poly
import matplotlib.pyplot as plt

ecutoff=800

x = np.linspace(0,ecutoff,50)
y = 1/np.sqrt(x+1)
y += 5*np.exp(-(x-0.7*ecutoff)**2/10)


plt.plot(x,y)

for n in [15, 30, 45, 60]:
    coef = poly.laguerre.lagfit(x,y,n)
    fit_y = poly.laguerre.lagval(x, coef)
    plt.plot(x, fit_y)
plt.show()
# xvals, weights = poly.leggauss(n)
