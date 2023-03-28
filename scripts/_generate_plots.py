import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter
import sys, traceback
import os.path as path
from QoL import set_highlighted_excepthook

set_highlighted_excepthook()


# Basic num arguments check
if  len(sys.argv) < 2:
    print("Usage: python generate_plots.py Carbon_1")
    print("Or to automatically use output mol file: python generate_plots.py Carbon_1 y")
    exit()
    

############
# File/directory names
#######  
dname_Figures = "../../../AC4DC_Figures/"
figures_ext = "" #.png
dname_Box_Sync = "~/Box Sync/" 
fname_charge_conservation = "charge_conservation"
fname_free = "free"
fname_HR_style = "HR_style"
fname_bound_dynamics = "bound_dynamics"  #Was called both _bound and dynamics so I assume it's bound dynamics or something.

# Plots to show
tot_charge = True
all_charge = True
free = True
#######

dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
print(dname_Figures)
label = sys.argv[1] +'_' 

#TODO make user inputs less big
if len(sys.argv) > 3:
     label += sys.argv[3] + '_'

name = sys.argv[1].replace('_',' ')
#name = ' '.join(re.split('(?<=.)(?=[A-Z])', sys.argv[1]))
#name = sys.argv[1].partition("_")[0] 

output_mol_query = ""
if len(sys.argv) > 2:
     output_mol_query = sys.argv[2]
pl = Plotter(sys.argv[1],output_mol_query)

if tot_charge: 
    pl.plot_tot_charge(every=10)
    pl.fig.set_size_inches(3,2.5)
    plt.savefig(dname_Figures + label + fname_charge_conservation + figures_ext)
    #plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_qcons.pgf')
    plt.close()

if free:  # Free graph breaks. 
    pl.fig.set_size_inches(3,2.5)
    pl.plot_free(log=True, min=1e-7, every=5)
    pl.fig_free.subplots_adjust(left=0.15)
    plt.savefig(dname_Figures + label + fname_free + figures_ext)
    plt.close()
    #plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_free.pgf')


if all_charge:
    pl.plot_all_charges()
    pl.fig.set_size_inches(3,2.5)
    pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
    plt.savefig(dname_Figures + label + fname_bound_dynamics + figures_ext)
    plt.close()

    #pl.print_bound_slice()
    #plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_bound.pgf')
    #plt.close()

    




#########
#  HR style
#####
#plt.close('all')



from plotter_core import fit_maxwell, maxwell

cmap = plt.get_cmap("tab10")

plt.rcParams["font.size"] = 16 # 8

''' #V2 DEFAULTS 
# slices = [-2.5, 0, 2.5, 5]
tt = 0
slices = [-7.5+tt, -5+tt, -2.5+tt, 0+tt]
# slices = [0, 15, 23]
energies = [200, 500, 500, 1000]  
'''

colrs = [cmap(i) for i in range(4)]

#Wiggle testing                       # Proportion of time steps (Starts at -16 fs unless otherwise stated)
#slices = [-15.5,-15,-14,-14]  #0.1
#slices = [-15.95,-15.85,-15.75,-15.65]  # 0.015 
#slices = [-29.90,-29.70,-29.50,-29.30]  # 0.015 15 fs
#slices = [-15.95,-15.95,-15.95,-15.95]  # 0.005
#slices = [-0.009,0,0.009,0.019,] #square keV 0.01 fs

#pulse wiggle testing
#slices = [-15.995,-15.98,-15.965,-15.95]  # 0.001 square
#slices = [-15.995,-15.98,-15.965,-15.95]  # 0.005
#slices = [-15.99,-15.9,-15.8,-15.8]  # 0.005 square
#slices = [-7.95,-7.95,-7.95,-7.95]  # square 

### shape, %, time,
#slices = [-9.5,-8.5,-7.5]      # square 10%, 10 fs,
slices = [-7.51]#[-9.8,-9.5] 
#slices = [-9.99,-9.98,-9.97]      # square << 1%, 10 fs, 

energies = [200,500,500,500,500,500]     # Cutoff energies for fitting MB curves. 
plot_fits = True # Whether to plot MB curves


for (t, e, col ) in zip(slices, energies, colrs):
        lines = pl.plot_step(t, normed=True, color = col, lw=0.5)
        if plot_fits:
            T = pl.plot_fit(t, e, normed=True, color=col, lw=0.5)  # Seems that e is the fitted cutoff energy for considering fitting. Default settings had it as approx. double the peak of the MB. -S.P.

# pl.fig_steps.set_size_inches(3,2.5)
pl.fig_steps.set_size_inches(6,5)

#all
pl.ax_steps.set_ylim([1e-4, 1])
#plt.yscale("linear")
#pl.ax_steps.set_ylim([-0.035, 0.035])
pl.ax_steps.set_xlim([1,10000]) 
#spike
#pl.ax_steps.set_ylim([1e-40, 1])
#pl.ax_steps.set_xlim([5000,7000]) 

pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

######### 
# Curve adding
#####
def add_curve(mat, fitE, **kwargs):
    fit = mat[:,0].searchsorted(fitE)
    Xdata = mat[:fit, 0]
    Ydata = mat[:fit, 1]/Xdata
    T, n = fit_maxwell(Xdata, Ydata)
    print(T, n)
    # Equivalent to plot_maxwell(T,n, **kwargs)
    X = np.logspace(0,4,100)
    pl.ax_steps.plot(X, maxwell(X, T, n)*X,
        '--', **kwargs)

# add_curve(X1, 300,  color = cmap(0),lw=0.5)
# add_curve(X2, 500,  color = cmap(1),lw=0.5)
# add_curve(X3, 500,  color = cmap(2),lw=0.5)
# add_curve(X4, 1000, color = cmap(3),lw=0.5)

#------\testing/-------
#pl.plot_maxwell(44.1,0.06*3/2)  # ~Sanders -7.5 fs - density of MB assumed to be 50% of total.  
pl.plot_maxwell(31,0.06*3/2,color = "r") #Hau-Riege
extr_handle,extr_label = pl.ax_steps.get_legend_handles_labels()
#########
# 
#####

handles, labels = pl.ax_steps.get_legend_handles_labels()
# e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  
order = list(range(0,len(labels) - 1,2)) + list(range(1,len(labels),2))

#if not one_fs_wiggle_testing:
pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left',ncol=2)
pl.ax_steps.legend(extr_handle,extr_label,ncol=2)

#plt.gcf()

plt.title(name + " - Free-electron distribution")
plt.tight_layout()

plt.savefig(dname_Figures + label + fname_HR_style + figures_ext)

print("Done! Remember to drink water!")

#TODO: Change size to match my screen by default, add --y option