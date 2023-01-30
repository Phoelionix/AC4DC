#%%
# fname = "Carbon_150_35"
# fname = "Sulfur_Sanders_Settings"
# fname = "C_180_60"


%matplotlib widget

#(use widget or ipympl)
import ipdb, pylab
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, RadioButtons
from plotter_core import Plotter #, fit_maxwell, maxwell
import pickle
import sys, traceback
from _save_interactive import WidgetHandler, set_highlighted_excepthook #from _save_interactive_copy import WidgetHandler_c, set_highlighted_excepthook_c
import os.path as path

#Freezes if use a function for some reason. Probably gets garbage collected or something.
'''
def run_interactive(interactive_directory):
    set_highlighted_excepthook()
    indir = interactive_directory # Definitely what indir means.
    #fname = sys.argv[1]
    #fname = input("Input file name: ")
    fname = "Carbon_150_35"
    fname = "Sulfur_Sanders_Settings"
    extension = ".fig.pickle"
    if len(fname) < len(extension) or fname[-len(extension):] != extension:
        fname += extension

    # Loads the pl object, 
    global pl
    pl = pickle.load(open(indir + fname, 'rb'))
    print(pl.fig_steps)
    plt.figure(pl.fig_steps)
    #%matplotlib widget

    widgies = WidgetHandler(pl)
    widget_update = widgies.update
    widget_update(0) 
    #plt.ion()
    plt.show()
indir = "../../../AC4DC_Interactives/"
indir = path.abspath(path.join(__file__ ,indir)) + "/"
run_interactive(indir)
'''

indir = "../../../AC4DC_Interactives/" # Definitely what indir means.
indir = path.abspath(path.join(__file__ ,indir)) + "/"
set_highlighted_excepthook()
#fname = sys.argv[1]
#fname = input("Input file name: ")
extension = ".fig.pickle"
if len(fname) < len(extension) or fname[-len(extension):] != extension:
    fname += extension

# Loads the pl object, 
#global pl
print("Loading pickle...")
pl = pickle.load(open(indir + fname, 'rb'))
print("Done!")
print(pl.fig_steps)
plt.figure(pl.fig_steps)
#%matplotlib widget

widgies = WidgetHandler(pl)
widget_update = widgies.update
widget_update(0) 
#plt.ion()
plt.show()


# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })


# def show_figure(fig):

#     # create a dummy figure and use its
#     # manager to display "fig"  
#     dummy = plt.figure()
#     new_manager = dummy.canvas.manager
#     new_manager.canvas.figure = fig
#     fig.set_canvas(new_manager.canvas)

# #show_figure(figx)

# #data = figx.axes[0].lines[0].get_data()
# # %%

# %%
fname = "C_180_60"
#fname = "C_215_100"
#fname = "C_215_50"
#fname = "C_250_100"
#fname = "C_215_80"
# %%
