import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D # <--- This is important for 3d plotting 
import sys

plt.rc('axes', labelsize=14)    # fontsize of the x and y labels

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
    fig, (ax1, ax2) = plt.subplots(nrows=1,ncols=2, sharex=True, sharey=True, figsize = (4,3))
    bottom = 0
    im = ax1.contourplot(energies, energies, -Qdata, vmin=bottom, vmax=bigg, cmap=cm)
    ax1.set_title(r'$-\sum_j Q^{TBR}_{j, kl}$', y=1.5)
    ax1.set_xlabel(r"$\epsilon_k$")
    ax1.set_ylabel(r"$\epsilon_l$")
    ax2.contoutplot(energies, energies, Gdata, vmin=bottom, vmax=bigg, cmap=cm)
    ax2.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl}$',y=1.5)
    ax2.set_xlabel(r"$\epsilon_k$")
    ax2.xaxis.set_label_coords(0.1, -0.05)
    ax1.xaxis.set_label_coords(0.1, -0.05)
    ax1.yaxis.set_label_coords(-0.05, 0.05)
    # ax3.pcolormesh(energies, energies, np.abs(Gdata+ Qdata), vmin=bottom, vmax=bigg, cmap=cm)
    # ax3.set_title(r'$\left|\sum_{\eta} \Gamma^{TBR}_{kl} + \sum_j Q^{TBR}_{j, kl}\right|$',y=1.05)
    # ax4.pcolormesh(energies, energies, Gdata+0.5*Qdata, vmin=bottom, vmax=bigg, cmap=cm)
    # ax4.set_title(r'$\sum_{\eta} \Gamma^{TBR}_{kl} + \frac{1}{2}\sum_j Q^{TBR}_{j, kl}$',y=1.05)
    
    fig.subplots_adjust(left=0.1, right=0.95,top=0.85, bottom=0.32)
    cbar_ax = fig.add_axes([0.15, 0.15, 0.7, 0.02]) 
    fig.colorbar(im, orientation="horizontal", cax=cbar_ax)
    


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
