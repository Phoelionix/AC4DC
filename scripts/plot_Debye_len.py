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

plt.rcParams["font.size"] = 8

slices = np.linspace(-15,15,100)

pl = Plotter(sys.argv[1])


tempList = []
denseList = []
for t in slices:
    tempList.append( pl.get_temp(t, 1000) ) # eV
    # denseList.append( pl.get_density(t) ) # per angstrom cube

T = np.array(tempList)
n = T [: ,1]
T = T[:,0]
print(n)

lambdaD=np.sqrt(C.epsilon_0 * C.angstrom * T *C.eV / n /C.e/C.e) # should have units Angstrom
slices = np.array(slices)
ax, ax2 = pl.setup_axes()
# ax.set_yscale('log')
ax.set_xlabel('Time (fs)')
ax.set_ylabel('Debye length (\AA)')


ax.plot(slices, lambdaD, label=r'$\lambda_D$')
ax2 = ax.twinx()
ax2.plot(slices, 4/3*C.pi*n*lambdaD**3, 'r--', label=r'$N_D$')
ax2.set_ylabel('Debye number $N_D$')
# ax.twinx().plot(slices, n, 'k:', label='$n$')
ax.legend()
pl.fig.subplots_adjust(left=0.15,bottom=0.2, top=0.95,right=0.85)
pl.fig.savefig('/Users/alaric-mba/Desktop/Debye_evolution.png')
# fig.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+sys.argv[1]+'_debye.pgf')
