import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import csv

raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)
plt.scatter(raw[:,1],raw[:,2],label=r'$\sum_\eta \Gamma^{EII}_{\xi \to \eta, k}$')
plt.scatter(raw[:,1],raw[:,3], label=r'$\sum_j \mathcal{Q}^{EII}_{\xi, jk}$')
plt.xlabel(r'$\bar{e}_k$')
plt.ylabel('Summed value (a.u.)')
plt.legend()
plt.show()
