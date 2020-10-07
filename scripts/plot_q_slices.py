import matplotlib
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter
import sys

slices = [-5, -2.5, 0, 2.5, 5]

pl = Plotter(sys.argv[1])
for t in slices:
    pl.plot_step(t)

plt.legend()
plt.savefig('/Users/alaric-mba/Desktop/free_distribution_evolution.png')
