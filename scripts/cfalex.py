import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter, get_colors
import sys
from labellines import *
import numpy as np

oldname='Carbon_HR2'
newname='Carbon_1um'

raw_old = np.genfromtxt('../code_from_alex/old__AC4DC/output/Charge_C.txt')
int_old = np.genfromtxt('../code_from_alex/old__AC4DC/output/Intensity_'+oldname+'.txt')

t=int_old[:,-1]
t2 = raw_old[:,-1]
old_charges =raw_old[:,:-1]
old_pulse = int_old[:,:-1]

fig,ax = plt.subplots(figsize=(3,2.5))
fig.subplots_adjust(bottom=0.2,right=0.95,top=0.95)
ax2 = ax.twinx()

x = 20

min_idx2 = t2.searchsorted(-x)
max_idx2 = t2.searchsorted(x)
min_idx = t.searchsorted(-x)
max_idx = t.searchsorted(x)



pl = Plotter(newname)
pl.aggregate_charges()
normed_charges = pl.chargeData['C'] / pl.chargeData['C'][0,0]


num_species = old_charges.shape[1]
# num_species = 5

for (i,c) in zip(range(num_species), get_colors(num_species,400)):
    ax.plot(pl.timeData, normed_charges[:,i], '-', linewidth=1, color=c)
    ax.plot(t2[min_idx2:max_idx2], old_charges[min_idx2:max_idx2,i], '--', label='%d+' % i, linewidth=1,color=c)
    
ax2.plot(t[min_idx:max_idx], old_pulse[min_idx:max_idx], 'k:', linewidth=1,alpha=0.5)
ax2.yaxis.set_visible(False)

ax.legend()
ax.set_xlabel('Time (fs)')

fig.savefig('/Users/alaric-mba/Desktop/cfalex.png')
outfile = '/Users/alaric-mba/Box Sync/Thesis/Figures/cfalex_'+sys.argv[1]+'.pgf'
fig.savefig(outfile)
print('saved figure to `' + outfile +'`')