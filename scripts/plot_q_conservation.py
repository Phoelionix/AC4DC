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

pl = Plotter(sys.argv[1])
# pl.plot_tot_charge()
# plt.savefig('/Users/alaric-mba/Desktop/charge_conservation.png')
pl.plot_free(max=1e-6)
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
