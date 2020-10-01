import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
import sys

rawQ = np.genfromtxt(sys.argv[1]+"_Q.csv")
rawGamma = np.genfromtxt(sys.argv[1]+"_gamma.csv")

energies = rawQ[:, 0]
Qdata = rawQ[:, 1:]
Gdata = rawGamma[:, 1:]


maxG = np.max(Gdata)
maxQ = np.max(Qdata)
bigg = max(maxQ,maxG)

#  colormap to use
cm = 'magma'

def plotit():    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2,ncols=2, sharex=True, sharey=True, figsize = (7,3))
    bottom = 0
    im = ax1.pcolormesh(energies, energies, -Qdata, vmin=bottom, vmax=bigg, cmap=cm)
    ax1.set_title(r'$-\sum_j Q^{TBR}_{j, kl}$', y=1.05)
    ax2.pcolormesh(energies, energies, Gdata, vmin=bottom, vmax=bigg, cmap=cm)
    ax2.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl}$',y=1.05)
    ax3.pcolormesh(energies, energies, np.abs(Gdata+ Qdata), vmin=bottom, vmax=bigg, cmap=cm)
    ax3.set_title(r'$\left|\sum_{\eta} \Gamma^{TBR}_{kl} + \sum_j Q^{TBR}_{j, kl}\right|$',y=1.05)
    ax4.pcolormesh(energies, energies, Gdata+0.5*Qdata, vmin=bottom, vmax=bigg, cmap=cm)
    ax4.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl} + \frac{1}{2}\sum_j Q^{TBR}_{j, kl}$',y=1.05)
    
    fig.subplots_adjust(right=0.8,top=0.85)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
    fig.colorbar(im, cax=cbar_ax, cmap=cm)
    


def plotsum():    
    fig, (ax) = plt.subplots(nrows=1,ncols=1, sharex=True, sharey=True, figsize = (3,3))
    im = ax.pcolormesh(energies, energies, (Gdata + Qdata), cmap = cm)
    ax.set_title('(G + Q)')
    fig.colorbar(im)

def plotdiff():    
    fig, (ax) = plt.subplots(nrows=1,ncols=1, sharex=True, sharey=True, figsize = (3,3))
    im = ax.pcolormesh(energies, energies, (Gdata + 0.5*Qdata), cmap = cm)
    ax.set_title('(G + 0.5Q)')
    fig.colorbar(im)

plotit()
# plotdiff()
# plotsum()
plt.show()
