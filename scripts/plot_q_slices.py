import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter
import sys

slices = [-5, -2.5, 0, 2.5]

pl = Plotter(sys.argv[1])
for t in slices:
    pl.plot_step(t, normed=False)

pl.ax_steps.set_ylim([1e-5, 5])
pl.ax_steps.set_xlim([1,10000])

plt.legend(['%1.1ffs' % i for i in slices])
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_6keV_no_ee.pgf')
