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
import numpy as np
import matplotlib.colors as colors

plee = Plotter('Carbon_square')
plne = Plotter('Carbon_square_no_ee')

fig, (ax1, ax2) = plt.subplots(1,2,sharey=True, sharex=True,figsize=(6.2,2.5))
fig.subplots_adjust(left=0.1, top=0.96, bottom=0.16,right=0.85)


min_n = 1e-7
max_n1 = np.max(plne.freeData.max())
max_n2 = np.max(plee.freeData.max())
max_n = max(max_n1, max_n2)

log = True

def plotit(ax, P):

    every = 1

    if every is None:
        Z = P.freeData.T
        T = P.timeData
    else:
        Z = P.freeData.T [:, ::every]
        T = P.timeData [::every]

    Z = np.ma.masked_where(Z < min_n, Z)
    Z = np.ma.masked_where(Z > max_n, Z)

    ax.set_facecolor('black')
    if log:
        cm = ax.pcolormesh(T, P.energyKnot*1e-3, Z, shading='gouraud',norm=colors.LogNorm(vmin=min_n, vmax=max_n),cmap='magma',rasterized=True)
    else:
        cm = ax.contourf(T, P.energyKnot, Z, 50, cmap='magma',rasterized=True)


    # plot the intensity
    ax_I = ax.twinx()
    ax_I.plot(T, P.intensityData[::every],'w:',linewidth=1)

    ax_I.get_yaxis().set_visible(False)
    return cm

plotit(ax1, plne)
cmap = plotit(ax2, plee)

ax1.set_ylabel("Energy (keV)")
ax1.set_xlabel("Time (fs)")
ax2.set_xlabel("Time (fs)")

cbar_ax = fig.add_axes([0.87, 0.16, 0.02, 0.80])

cbar = fig.colorbar(cmap, cax=cbar_ax)
cbar.ax.set_ylabel('Free Electron Density, Å$^{-3}$', rotation=270,labelpad=10)

fig.savefig('/Users/alaric-mba/Desktop/eefree.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/compare_CHR2_free.pgf')