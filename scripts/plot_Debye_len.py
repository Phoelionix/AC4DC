import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter, fit_maxwell
import sys
import numpy as np
from scipy import constants as C

# from labellines import *

plt.rcParams["font.size"] = 16

pl = Plotter(sys.argv[1])
slices = np.linspace(pl.timeData[0],pl.timeData[-1],50)

tempList = []
denseList = []
for t in slices:
    tempList.append( pl.get_temp(t, 1000) ) # eV
    # denseList.append( pl.get_density(t) ) # per angstrom cube

T = np.array(tempList)
n = T [: ,1]
T = T[:,0]
print(n)
print(T)

lambdaD=np.sqrt(C.epsilon_0 * C.angstrom * T *C.eV / n /C.e/C.e) # should have units Angstrom
slices = np.array(slices)
fig, (ax, ax2) = plt.subplots(nrows=2, sharex=True, figsize=(6,5))
# ax.set_yscale('log')
ax.set_ylabel(r'Debye length $\lambda_D$ (\AA)')
ax.set_ylim([0,5])
ax.plot(slices, lambdaD, label=r'$\lambda_D$')
ax_int = ax.twinx()
ax_int.plot(pl.timeData, pl.intensityData, lw = 1, c = 'black', ls = ':', alpha = 0.7)
# ax_int.set_ylabel('Pulse Intensity (photons cm$^{-2}$ s$^{-1}$)')
ax2.set_xlabel("Time (fs)")
ax_int.axes.get_yaxis().set_visible(False)

ax2.set_yscale('log')
ax2.get_yaxis().get_major_formatter().labelOnlyBase = False
ax2.plot(slices, 4/3*C.pi*n*lambdaD**3, 'r--', label=r'$N_D$')
ax2.set_ylabel('Debye number $N_D$')
ax2.set_ylim([0.1,1e4])

ax_int = ax2.twinx()
ax_int.plot(pl.timeData, pl.intensityData, lw = 1, c = 'black', ls = ':', alpha = 0.7)
# ax_int.set_ylabel('Pulse Intensity (photons cm$^{-2}$ s$^{-1}$)')
ax2.set_xlabel("Time (fs)")
ax_int.axes.get_yaxis().set_visible(False)

fig.subplots_adjust(left=0.2,bottom=0.15, top=0.95,right=0.9)
fig.savefig('/Users/alaric-mba/Desktop/Debye_evolution.png')
fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_debye.pgf')
