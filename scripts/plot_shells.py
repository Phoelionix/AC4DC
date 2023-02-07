#%%
# For comparison with Hau-Riege 2007. 10.1103/PhysRevA.76.042511
%matplotlib widget

import matplotlib
# matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter, parse_elecs_from_latex
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np



#handle = "Carbon_Sanders"                                       
#atoms = [('C' , ['1s', '2s','2p'], [2,2,2])]                    

handle = "Carbon_Wiggleless" 
atoms = [('C_fast' , ['1s', '2N'], [2,4])]





plt.rcParams["font.size"] = 8


label = handle.replace('_',' ')
pl = Plotter(handle,"y");
# close step axes
plt.close()


def get_tot_occ(a, subshell, start_charge = None, n=None):

    states = pl.statedict[a]

    norm = [1 for elem in states]
    if start_charge != None:
        # Y axis becomes average charge
        ax.set_ylabel(r"$\overline{Z}$",rotation="horizontal",size=12)
        ax.yaxis.set_label_coords(-.1, .5)
        orboccs = parse_elecs_from_latex(states[0])
        norm = start_charge/(orboccs[subshell]*pl.boundData[a][0,0])

    tot = np.zeros_like(pl.boundData[a][:,0])
    for i in range(len(states)):
        orboccs = parse_elecs_from_latex(states[i])
        if n is None:
            tot += norm*orboccs[subshell]*pl.boundData[a][:,i]
        elif orboccs[subshell] == n:
             tot += norm*pl.boundData[a][:,i]

    return tot

# fig, ax = plt.subplots(figsize=(3,2.5))
ax, _x = pl.setup_axes()
# ax.set_title("Charge state dynamics")

ax.set_ylabel(r"Electron density (\AA$^{-3}$)")

pl.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)

# inset axes....
# axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

y=[]
for atom, shells, charges in atoms: 
    for i in range(len(shells)):
        shell, charge = shells[i], charges[i]
        shell_substitute = shell
        # Shell
        if shell[-1] == 'N': 
            shell_substitute = shell[:-1]+ 'p'  # Orbitals are approximated as p orbitals after slater energy. At least as of coding.
            y = get_tot_occ(atom, shell_substitute, start_charge = charge)
            label = atom + ' ' + shell
        # Subshells, go through each.
        else:
            subsh_occ = get_tot_occ(atom, shell, start_charge = charge)
            if len(y) == 0:
                y = subsh_occ
                label = atom + ' ' + shell
            else:
                y += subsh_occ
                label += ' ' + shell
            if i != len(shells) - 1 and shells[i+1][-1] != 's' and shells[i+1][-1] != 'N':
                continue
        ax.plot(pl.timeData, y, label= label) 
        y=[]
# axbox = ax.get_position()

# pl.fig.legend( loc = (axbox.x0 + 0.5, axbox.y0 + 0.15))
ax.legend()

ax.xaxis.set_minor_locator(MultipleLocator(2))
ax.yaxis.set_major_locator(MultipleLocator(1))  #ax.yaxis.set_minor_locator(MultipleLocator(0.004))
ax.set_ylim(0,4.04)
ax.set_xlim(-10,-7.5)
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
# %%
