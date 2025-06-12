#%%
%matplotlib widget


import os
os.getcwd() 
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
print(sys.path)


import matplotlib
#matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 8,
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
from plotter_core import Plotter
import sys, traceback
import os.path as path
from QoL import set_highlighted_excepthook
import numpy as np


set_highlighted_excepthook()

DPI=800
FIGWIDTH = 3.34646
FIGHEIGHT = FIGWIDTH*2/3
###
ANG_TO_BOHR =0.529177
HIGH_RES = 1
LOW_RES = 4 # Used only if RESOLUTION=TRUE
NUM_TRACES = 5
SHOW_AVG=False
SHOW_SNAPSHOTS=True
PLOT_DIFF = True
RESOLUTION= True
YLIM = [-35,5]
#YLIM = [1.5,6]
###
#handles = ["lys_nass_no_S_3","lys_nass_gauss","lys_nass_Gd_full_1"]
#handles=["lys_nass_no_S_3"]
#handles = ["I3C_55fs_4"]
#handles = ["carbon_ff_test_15_3"]; start_t = -18; end_t = 18
#handles = ["carbon_ff_test_3_2"]; start_t = -3*1.2; end_t = 3*1.2
handles = ["CB_11","lys_Gd_salt_heavy_atoms_off_4","lys_solvated_light_10"]; start_t = -18; end_t = 17.98
##
#atoms = pl.atomdict
#atoms = ["N","Gd_fast"]
#atoms = ["C","N","O"]
#atoms = ["I_fast"]
atoms = ["C"]
#######
dname_Figures = "../../output/_Graphs/plots/"
dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
for handle in handles:
    label = handle +'_' 

    name = handle.replace('_',' ')

    pl = Plotter(handle)
    plt.close()
    photon_energy = 6000
    # q_max is in units of bohr^-1 (atomic units). 
    q_max = 2*np.pi/HIGH_RES*ANG_TO_BOHR
    pl.initialise_form_factor_params(start_t,end_t,q_max,photon_energy,50,100)


    if RESOLUTION:
        q_min = 2*np.pi/LOW_RES*ANG_TO_BOHR
    else:
        q_min=0
    for atom in atoms:
        # Plot atomic form factor at different times as a function of q
        #pl.plot_ffactor_get_R_AC4DC(atom)   # More accurate form factors from AC4DC

        # Slater form factors. Should be fine for light atoms.

        pl.plot_form_factor(NUM_TRACES,[atom],resolution=RESOLUTION,q_min=q_min,fig_width=FIGWIDTH,fig_height=FIGHEIGHT,show_average=SHOW_AVG,show_snapshots=SHOW_SNAPSHOTS,percentage_change=PLOT_DIFF)
        #plt.title(handle+" - "+atom) 
        plt.tight_layout()
        plt.ylim(YLIM)
        plt.savefig(dname_Figures+atom+"-ff-"+handle+".png",dpi=DPI)
    if len(atoms) > 1:
        pl.plot_form_factor(NUM_TRACES,atoms,resolution=RESOLUTION,q_min=2*np.pi/LOW_RES*ANG_TO_BOHR,fig_width=FIGWIDTH,fig_height=FIGHEIGHT,show_average=SHOW_AVG,show_snapshots=SHOW_SNAPSHOTS,percentage_change=PLOT_DIFF)
        plt.title(handle+" - combined") 

    #print(pl.get_A_bar(-10,-7.5,12,12,"C_fast","C_fast",100))

    vmin = 0
    vmax = 1

    #pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax)
    #pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax,title = r"Naive Foreground coherence $\bar{A}$")
    #pl.plot_R_factor()
    #pl.print_bound_slice(-10)
    #print(" GAP ")
    #pl.print_bound_slice(-7.5)

    print("Done! Remember to stretch!")

# %%
