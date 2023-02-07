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

handle =  "Carbon_Wiggleless"#"Carbon_Sanders"
#######

label = handle +'_' 

name = handle.replace('_',' ')

pl = Plotter(handle,"y")
plt.close()

q_max = 2*np.pi#2*2*np.pi for sanders  # HR significant range is 0.1 to 0.6 lattice spacings 2*np.pi*0.6
# form factors dependent on q 
#pl.plot_ffactor_time_slices("C_fast")
# Atomic form factor at time as a function of q

pl.plot_form_factor(np.linspace(-10,-7.5,11),q_max) 

print(pl.get_A_bar(-10,-7.5,12,12,"C_fast","C_fast",100))

vmin = 0
vmax = 1
pl.plot_A_map(-10,-7.5,q_max,"C_fast","C_fast",100,100,vmin=vmin,vmax=vmax)
pl.plot_A_map(-10,-7.5,q_max,"C_fast","C_fast",100,100,vmin=vmin,vmax=vmax,naive=True,title = r"Naive Foreground coherence $\bar{A}$")

#pl.print_bound_slice(-10)
#print(" GAP ")
#pl.print_bound_slice(-7.5)

print("Done! Remember to stretch!")

# %%
