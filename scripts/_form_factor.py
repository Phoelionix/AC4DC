#%%
%matplotlib widget
import matplotlib
#matplotlib.use("pgf")
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
import numpy as np


set_highlighted_excepthook()

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

handle =  "Carbon_294"#"Carbon_Wiggleless_2500_Fluence"#"Carbon_Sanders"""
#######

label = handle +'_' 

name = handle.replace('_',' ')

pl = Plotter(handle,"y")
plt.close()
photon_energy = 6000
q_max = 2#pl.theta_to_q(22,photon_energy) # In units of bohr^-1. 
pl.initialise_coherence_params(-10,-6.4,q_max,photon_energy,20,100,True)




# form factors dependent on q 
pl.plot_ffactor_time_slices("C_fast")
# Atomic form factor at time as a function of q



pl.plot_form_factor(np.linspace(-10,-6.4,11)) 

#print(pl.get_A_bar(-10,-7.5,12,12,"C_fast","C_fast",100))

vmin = 0
vmax = 1

#pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax)
#pl.plot_A_map("C_fast","C_fast",vmin=vmin,vmax=vmax,title = r"Naive Foreground coherence $\bar{A}$")
pl.plot_R_factor()
#pl.print_bound_slice(-10)
#print(" GAP ")
#pl.print_bound_slice(-7.5)

print("Done! Remember to stretch!")

# %%
