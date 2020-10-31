import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import fit_maxwell, maxwell

plt.rcParams["font.size"] = 8
cmap = plt.get_cmap("tab10")

import sys
import numpy as np

raw = np.genfromtxt(sys.argv[1],skip_header=2,delimiter=',')

fig,ax = plt.subplots()

fig.set_size_inches(3,2.5)
# fig.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)
fig.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)

X1 = raw[:,0:2]
X2 = raw[:,2:4]
X3 = raw[:,4:6]
X4 = raw[:,6:8]

X1 = X1[X1[:,0].argsort()]
X2 = X2[X2[:,0].argsort()]
X3 = X3[X3[:,0].argsort()]
X4 = X4[X4[:,0].argsort()]

ax.plot(X1[:,0],X1[:,1],'.', label='-7.5fs, T=31eV',  color = cmap(0), markersize=0.5)
ax.plot(X2[:,0],X2[:,1],'.', label='-5fs, T=70eV',    color = cmap(1), markersize=0.5)
ax.plot(X3[:,0],X3[:,1],'.', label='-2.5fs, T=125eV', color = cmap(2), markersize=0.5)
ax.plot(X4[:,0],X4[:,1],'.', label='0fs, T=195eV',    color = cmap(3), markersize=0.5)

def add_curve(mat, fitE, **kwargs):
    fit = mat[:,0].searchsorted(fitE)
    Xdata = mat[:fit, 0]
    Ydata = mat[:fit, 1]/Xdata
    T, n = fit_maxwell(Xdata, Ydata)
    print(T, n)
    X = np.logspace(0,4,100)
    ax.plot(X, maxwell(X, T, n)*X,
        '--', **kwargs)

add_curve(X1, 300,  color = cmap(0),lw=0.5)
add_curve(X2, 500,  color = cmap(1),lw=0.5)
add_curve(X3, 500,  color = cmap(2),lw=0.5)
add_curve(X4, 1000, color = cmap(3),lw=0.5)

ax.set_ylim([1e-4, 1]) 
ax.set_xlim([1,10000]) 
ax.axhline(y=1e-3,color='k', linestyle='--', lw=0.5)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel(r'$f(\epsilon) \Delta \epsilon$')

ax.legend()
plt.savefig('/Users/alaric-mba/Desktop/HR_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/HauRiege_slices.pgf')