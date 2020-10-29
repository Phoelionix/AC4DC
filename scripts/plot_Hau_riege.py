import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 8

import sys
import numpy as np

raw = np.genfromtxt(sys.argv[1],skip_header=2,delimiter=',')

fig,ax = plt.subplots()

fig.set_size_inches(3,2.5)
fig.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)

X1 = raw[:,0:2]
X2 = raw[:,2:4]
X3 = raw[:,4:6]
X4 = raw[:,6:8]

X1 = X1[X1[:,0].argsort()]
X2 = X2[X2[:,0].argsort()]
X3 = X3[X3[:,0].argsort()]
X4 = X4[X4[:,0].argsort()]

ax.plot(X1[:,0],X1[:,1], label='-7.5fs, T=31eV',lw=0.5)
ax.plot(X2[:,0],X2[:,1], label='-5fs, T=70eV',lw=0.5)
ax.plot(X3[:,0],X3[:,1], label='-2.5fs, T=125eV',lw=0.5)
ax.plot(X4[:,0],X4[:,1], label='0fs, T=195eV',lw=0.5)

ax.set_ylim([1e-3, 2e-1]) 
ax.set_xlim([1,10000]) 

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy, eV')
ax.set_ylabel(r'$f(\epsilon) \Delta \epsilon$')

ax.legend()
plt.savefig('/Users/alaric-mba/Desktop/HR_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/HauRiege_slices.pgf')