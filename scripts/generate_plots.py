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

plt.rcParams["font.size"] = 9

slices = [-5, -2.5, 0, 2.5]
colours = ['red', 'green', 'blue', 'purple']

pl = Plotter(sys.argv[1])
for (t, c) in zip(slices, colours):
    pl.plot_step(t, normed=True, fitE=2000)

pl.ax_steps.set_ylim([1e-5, 1])
pl.ax_steps.set_xlim([1,10000])
pl.fig_steps.subplots_adjust(bottom=0.3)
pl.fig_steps.legend(loc='lower left', ncol=4)
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_6keV_no_ee.pgf')

pl.plot_tot_charge()
plt.savefig('/Users/alaric-mba/Desktop/charge_conservation.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_qcons_6keV_no_ee.pgf')

pl.plot_free(log=True, )
plt.savefig('/Users/alaric-mba/Desktop/charge_conservation.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_free_6keV_no_ee.pgf')
