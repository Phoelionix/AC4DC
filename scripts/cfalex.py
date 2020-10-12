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

raw_old = np.genfromtxt('../code_from_alex/old__AC4DC/output/Charge_C.txt')
int_old = np.genfromtxt('../code_from_alex/old__AC4DC/output/Intensity_Carbon_HR2.txt')

t=int_old[:,-1]
t2 = raw_old[:,-1]
old_charges =raw_old[:,:-1]
old_pulse = int_old[:,:-1]

fig,ax = plt.subplots(figsize=(3,2.5))
fig.subplots_adjust(bottom=0.2,right=0.95,top=0.95)
ax2 = ax.twinx()
min_idx = t.searchsorted(-20)
max_idx = t.searchsorted(20)



pl = Plotter('Carbon_HR2')
pl.aggregate_charges()
normed_charges = pl.chargeData['C'] / pl.chargeData['C'][0,0]


num_species = old_charges.shape[1]

for (i,c) in zip(range(num_species), get_colors(num_species,400)):
    l = ax.plot(pl.timeData, normed_charges, '-', linewidth=1, alpha=0.8)
    ax.plot(t[min_idx:max_idx], old_charges[min_idx:max_idx,i], '--', label='%d+' % i, linewidth=1,alpha=0.8,color=l[-1].get_color())
    
ax2.plot(t2[min_idx:max_idx], old_pulse[min_idx:max_idx], ':', linewidth=1,alpha=0.5)
ax2.yaxis.set_visible(False)

ax.legend()
ax.set_xlabel('Time (fs)')

fig.savefig('/Users/alaric-mba/Desktop/cfalex.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/cfalex_6keV_1000nm.pgf')