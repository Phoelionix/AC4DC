import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import csv

raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)
fig, ax = plt.subplots()
ax.scatter(raw[:,1],raw[:,2],label=r'$\sum_\eta \Gamma^{EII}_{\xi \to \eta, k}$')
ax.scatter(raw[:,1],raw[:,3], label=r'$\sum_j \mathcal{Q}^{EII}_{\xi, jk}$')
ax.set_xlabel(r'Mean energy at index $k$ (Ha)')
ax.set_ylabel('Summed value (a.u.)')

ax2 = ax.twinx()

ax2.set_ylabel('Absolute difference')
ax2.plot(raw[:,1],np.abs(raw[:,2]-raw[:,3]),'k.',label='Absolute Difference')
ax2.set_yscale('log')

ax.legend(loc=0)
ax2.legend(loc='lower right')
ax2.grid()
plt.show()
