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
from matplotlib.ticker import MultipleLocator
import numpy as np

plt.rcParams["font.size"] = 8

handle = 'water_strong'

label = handle

pl = Plotter(handle)

def get_tot_occ(a, subshell, n=None):
    states = pl.statedict[a]
    tot = np.zeros_like(pl.boundData[a][:,0])
    for i in range(len(states)):
        orboccs = parse_elecs_from_latex(states[i])
        if n is None:
            tot += orboccs[subshell]*pl.boundData[a][:,i]
        elif orboccs[subshell] == n:
             tot += pl.boundData[a][:,i]

    return tot

# fig, ax = plt.subplots(figsize=(3,2.5))
ax, _x = pl.setup_axes()

# ax.set_title("Charge state dynamics")
ax.set_ylabel(r"Electron density (\AA$^{-3}$)")

pl.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)

# inset axes....
# axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

a = 'O'

for subsh in ['1s', '2s', '2p']:
    # for n in range(3):
        # name = '$%s^%d$' % (subsh, n)
    ax.plot(pl.timeData, get_tot_occ(a,subsh), label='O ' + subsh)
    # axins.plot(pl.timeData, get_tot_occ(a,subsh), label='O' + subsh)

ax.plot(pl.timeData, get_tot_occ('H','1s'), label='H 1s')

# axbox = ax.get_position()

# pl.fig.legend( loc = (axbox.x0 + 0.5, axbox.y0 + 0.15))
ax.legend()

ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.yaxis.set_minor_locator(MultipleLocator(0.004))
ax.tick_params(which='both', direction='in')

# # sub region of the original image
# x1, x2, y1, y2 = -10.1, -8, 0.093 , 0.112
# axins.set_xlim(x1, x2)
# axins.set_ylim(y1, y2)
# axins.set_xticklabels('')
# axins.set_yticklabels('')

# ax.legend()
# ax.set_xlabel('Time (fs)')


# ax.indicate_inset_zoom(axins)

plt.savefig('/Users/alaric-mba/Desktop/subshell_evolution.png')
plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_subshells.pgf')