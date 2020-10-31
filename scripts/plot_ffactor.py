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

label = sys.argv[1] +'_' + sys.argv[2]

pl = Plotter(sys.argv[1])

fig, ax = pl.plot_ffactor('C', 10, timespan=(-10,0))
sm = plt.cm.ScalarMappable(cmap=plt.get_cmap('plasma'), norm=plt.Normalize(vmin=-10, vmax=0))
 
cbaxes = fig.add_axes([0.5, 0.9, 0.35, 0.02])
a = plt.colorbar(sm, orientation='horizontal', cax = cbaxes)
cbaxes.set_xlabel('Time (fs)')

fig.set_size_inches(3,2.5)
fig.subplots_adjust(left=0.2, bottom=0.2 ,right=0.95, top=0.95)
fig.savefig('/Users/alaric-mba/Desktop/ffactor.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_ffactor.pgf')
