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
    'font.size': 10,
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

#handle =  "lys_nass_no_S_2"#"sulfur_3shell_baseline_3"#"Naive_Lys_mid_21"#"Carbon_294"#"Carbon_Wiggleless_2500_Fluence"#"Carbon_Sanders"""
#handle =  "lys_nass_6"#
#handle =  "lys_nass_Gd_12"#
handles = ["lys_nass_no_S_2","lys_nass_6","lys_nass_Gd_16"]
#######
for handle in handles:
    label = handle +'_' 

    name = handle.replace('_',' ')

    pl = Plotter(handle)
    plt.close()
    photon_energy = 6000
    # q_max is in units of bohr^-1 (atomic units). 
    q_max = 20*0.529177# 20 angstrom. 
    pl.initialise_form_factor_params(-10,17.5,q_max,photon_energy,50,100,True)

    # form factors dependent on q 
    #atoms = pl.atomdict
    atoms = ["N"]
    for atom in atoms:
        pl.plot_ffactor_get_R_sanders(atom)
        # Atomic form factor at time as a function of q



        pl.plot_form_factor(6,[atom])
        plt.title(handle+" - "+atom) 
    #pl.plot_form_factor(6)
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
