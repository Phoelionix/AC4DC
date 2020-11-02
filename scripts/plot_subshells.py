import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter, parse_elecs_from_latex
import sys
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["font.size"] = 8

handle = 'Carbon_HR1'

label = handle

pl = Plotter(handle)

def get_tot_occ(a, subshell):
    states = pl.statedict[a]
    tot = np.zeros_like(pl.boundData[a][:,0])
    for i in range(len(states)):
        orboccs = parse_elecs_from_latex(states[i])
        tot += orboccs[subshell]*pl.boundData[a][:,i]
    return tot

# fig, ax = plt.subplots(figsize=(3,2.5))
ax, _x = pl.setup_axes()

# ax.set_title("Charge state dynamics")
ax.set_ylabel(r"Electron density (\AA$^{-3}$)")

pl.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)

# inset axes....
axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

def plot_subshell(a, subsh):
    ax.plot(pl.timeData, get_tot_occ(a,subsh), label=subsh)

def plot_subshell_zoomed(a, subsh):
    axins.plot(pl.timeData, get_tot_occ(a,subsh), label=subsh)

for a, s in zip(['C', 'C', 'C'], ['1s', '2s', '2p']):
    plot_subshell( a, s)
    plot_subshell_zoomed( a, s)


axbox = ax.get_position()

pl.fig.legend( loc = (axbox.x0 + 0.5, axbox.y0 + 0.15))

# sub region of the original image
x1, x2, y1, y2 = -10.1, -8, 0.1, 0.113
axins.set_xlim(x1, x2)
axins.set_ylim(y1, y2)
axins.set_xticklabels('')
axins.set_yticklabels('')

# ax.legend()
# ax.set_xlabel('Time (fs)')


ax.indicate_inset_zoom(axins)

plt.savefig('/Users/alaric-mba/Desktop/subshell_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_subshells.pgf')