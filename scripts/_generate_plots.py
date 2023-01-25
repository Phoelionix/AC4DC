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
import sys

############
# Basic num arguments check
#######  
# https://stackoverflow.com/questions/14775916/coloring-exceptions-from-python-on-a-terminal
def set_highlighted_excepthook():  
    import sys, traceback
    from pygments import highlight
    from pygments.lexers import get_lexer_by_name
    from pygments.formatters import TerminalFormatter

    lexer = get_lexer_by_name("pytb" if sys.version_info.major < 3 else "py3tb")
    formatter = TerminalFormatter()

    def myexcepthook(type, value, tb):
        tbtext = ''.join(traceback.format_exception(type, value, tb))
        sys.stderr.write(highlight(tbtext, lexer, formatter))

    sys.excepthook = myexcepthook

set_highlighted_excepthook()
if  len(sys.argv) != 3:
    raise Exception("Usage: python ./_generate_plots.py Carbon_HR1 NAME")

############
# File/directory names
#######  
dname_Figures = "../../AC4DC_Figures/"
figures_ext = "" #.png
dname_Box_Sync = "~/Box Sync/" 
fname_charge_conservation = "charge_conservation"
fname_free = "free"
fname_HR_style = "HR_style"
fname_bound_dynamics = "bound_dynamics"  #Was called both _bound and dynamics so I assume it's bound dynamics or something.

one_fs_wiggle_testing = True
#######


label = sys.argv[1] +'_' + sys.argv[2] + '_'

pl = Plotter(sys.argv[1])

#if not one_fs_wiggle_testing: # don't need these plots now
pl.plot_tot_charge(every=10)
pl.fig.set_size_inches(3,2.5)
plt.savefig(dname_Figures + label + fname_charge_conservation + figures_ext)
#plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_qcons.pgf')


if not one_fs_wiggle_testing:  # Free graph breaks. 
    pl.fig.set_size_inches(3,2.5)
    pl.plot_free(log=True, min=1e-7, every=5)
    pl.fig_free.subplots_adjust(left=0.15)
    plt.savefig(dname_Figures + label + fname_free + figures_ext)
    plt.close()
    #plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_free.pgf')



pl.plot_all_charges()
pl.fig.set_size_inches(3,2.5)
pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
plt.savefig(dname_Figures + label + fname_bound_dynamics + figures_ext)
#plt.savefig('/Users/alaric-mba/Box Sync/Thesis/Figures/'+label+'_bound.pgf')

plt.close()
plt.close()

#########
#  HR style
#####
#plt.close('all')



from plotter_core import fit_maxwell, maxwell

cmap = plt.get_cmap("tab10")

plt.rcParams["font.size"] = 16 # 8

# slices = [-2.5, 0, 2.5, 5]
tt = 0
slices = [-7.5+tt, -5+tt, -2.5+tt, 0+tt]
# slices = [0, 15, 23]
energies = [200, 500, 500, 1000]  #MB upper boundary on energy? -S.P.
colrs = [cmap(i) for i in range(4)]

#Wiggle testing                       # Proportion of time steps (Starts at -16 fs unless otherwise stated)
#wiggle_slices = [-15.5,-15,-14,-14]  #0.1
#wiggle_slices = [-15.95,-15.85,-15.75,-15.65]  # 0.015 
#wiggle_slices = [-29.90,-29.70,-29.50,-29.30]  # 0.015 15 fs
#wiggle_slices = [-15.95,-15.95,-15.95,-15.95]  # 0.005
#wiggle_slices = [-0.009,0,0.009,0.019,] #square keV 0.01 fs

#pulse wiggle testing
#wiggle_slices = [-15.995,-15.98,-15.965,-15.95]  # 0.001 square
#wiggle_slices = [-15.995,-15.98,-15.965,-15.95]  # 0.005
#wiggle_slices = [-15.99,-15.9,-15.8,-15.8]  # 0.005 square
#wiggle_slices = [-7.95,-7.95,-7.95,-7.95]  # square 

### shape, %, time,
wiggle_slices = [-9.5,-8.5,-7.5]      # square 10%, 10 fs, 
#wiggle_slices = [-9.99,-9.98,-9.97]      # square << 1%, 10 fs, 

wiggle_energies = [500,500,500]
#wiggle_energies = [1000]
#wiggle_slices = [-7.5]

if not one_fs_wiggle_testing:  
    for (t, e, col ) in zip(slices, energies, colrs):
        lines = pl.plot_step(t, normed=True, color = col, lw=0.5)
        T = pl.plot_fit(t, e, normed=True, color=col, lw=0.5)
else:
    for (t, e, col ) in zip(wiggle_slices, wiggle_energies, colrs):
        lines = pl.plot_step(t, normed=True, color = col, lw=0.5)
        T = pl.plot_fit(t, e, normed=True, color=col, lw=0.5)  # Seems that e is the cutoff energy. Approx. double the peak of the MB.

# pl.fig_steps.set_size_inches(3,2.5)
pl.fig_steps.set_size_inches(6,5)

#all
pl.ax_steps.set_ylim([1e-4, 1])
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
#########
# 
#####

handles, labels = pl.ax_steps.get_legend_handles_labels()
# e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  
order = list(range(0,len(labels) - 1,2)) + list(range(1,len(labels),2))

#if not one_fs_wiggle_testing:
pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper left',ncol=2)
#plt.gcf()
plt.savefig(dname_Figures + label + fname_HR_style + figures_ext)

print("Done! Remember to drink water!")