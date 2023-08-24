import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter
import sys
import numpy as np


runs_to_compare = ['Carbon_HR1', 'Carbon_HR3', 'Carbon_HR2']
names = ['100nm', '500nm', '1000nm']

plotters = []

for f in runs_to_compare:
    p = Plotter(f)
    p.aggregate_charges()
    plotters.append(p)


every=10

ax, ax2 = plotters[0].setup_axes()
plotters[0].fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)
for p, name in zip(plotters, names):
    T = p.timeData[::every]
    Q = np.zeros(T.shape[0])
    for a in p.atomdict:
        atomic_charge = np.zeros(T.shape[0])
        for i in range(p.chargeData[a].shape[1]):
            atomic_charge += p.chargeData[a][::every,i]*i
        # ax.plot(T, atomic_charge, label = a)
        Q += atomic_charge

    de = np.append(p.energyKnot, p.energyKnot[-1]*2 - p.energyKnot[-2])
    de = de [1:] - de[:-1]
    tot_free_Q =np.dot(p.freeData, de)

    # Q now stores total (+ve) ionic charge, tot_free_Q stores free charge
    lines = ax.plot(T, Q, '-', label=name)
    ax.plot(T,tot_free_Q[::every] , ':', color=lines[-1].get_color()) 
    # plt.plot(T,tot_free_Q, ':')

ax.legend()
ax.set_ylabel(r'Density, \AA$^{-3}$')

plt.savefig('/Users/alaric-mba/Desktop/charge_conservation.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/Carbon_HR12_qcons_6keV_no_ee.pgf')
