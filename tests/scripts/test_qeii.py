import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

plt.style.use('seaborn-muted')

import sys
import csv


raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)
fig, ax = plt.subplots(figsize = (4,3))
plt.subplots_adjust(left=0.16,bottom=0.16,right=0.95,top=0.95)
y1 = raw[:,2]
y2 = raw[:,3]
ax.plot(raw[:,1],np.ma.masked_where(y1<=0, y1),'.',label=r'$\sum_\eta \Gamma^{EII}_{\xi \to \eta, k}$')
ax.plot(raw[:,1],0.1*np.ma.masked_where(y2<=0, y2),'.', label=r'$0.1 \times \sum_j \mathcal{Q}^{EII}_{\xi, jk}$')
ax.set_xlabel(r'Mean energy at index $k$ (Ha)')
ax.set_ylabel(r'$R_{\xi, k}$')
Y = np.abs(raw[:,2]-raw[:,3])
Y = np.ma.masked_where(Y<=0,Y)

ax.tick_params(axis='y',which='minor')
ax.set_ylim([1e-6,1e2])
locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.2,0.4,0.6,0.8),numticks=12)
ax.yaxis.set_minor_locator(locmin)
ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.plot(raw[:,1],Y,'k.',label='Residuals')


ax.legend(loc=0)
ax.set_yscale('log')

fig.savefig('/Users/alaric-mba/Desktop/eii.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/Qeii_conservation.pgf')
