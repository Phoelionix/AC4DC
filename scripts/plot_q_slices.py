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
from labellines import *

plt.rcParams["font.size"] = 8

slices = [-5, -2.5, 0, 2.5]
colours = ['red', 'green', 'blue', 'purple']

pl = Plotter(sys.argv[1])

lines = []
fitlines = []


for (t, c) in zip(slices, colours):
    lines.append(pl.plot_step(t, normed=True)[-1])
    fitlines.append(pl.plot_fit(t, 2000, normed=True, c=lines[-1].get_color())[-1])

pl.fig_steps.set_size_inches(3,2.5)

pl.ax_steps.set_ylim([1e-5, 10])
pl.ax_steps.set_xlim([1,10000])

pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.85,top=0.95)
pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

pl.fig_steps.legend(loc='lower center', ncol=4,labelspacing=0)
# labelLines(fitlines,align=False)
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_6keV_no_ee.pgf')
