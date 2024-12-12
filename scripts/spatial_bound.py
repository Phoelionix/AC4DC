import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'font.size': 10,
    'text.usetex': True,
    'pgf.rcfonts': False,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook

CHARGE_DIFFERENCE = False # Set True if want ionised atoms' trace to start from origin
PLOT_DERIVATIVE = False # Plot the rate of avg charge gain
PLOT_MODE = 1  # 0: plot all charges, 1: plot element total charges # 2: plot orbital charges
YLIM = [0,None]
FIGWIDTH = 3.49751*2/3
#FIGHEIGHT = 3.49751/2
PLOTHEIGHT = 3.49751/2
#XLIM = [None,None]
XLIM = [None,None]
#YLIM=[0,6]

CUSTOM_LEGEND=("1s: CNO","2s: CNO","2p: CNO","1s: CNO,S,Gd","2s: CNO,S,Gd","2p: CNO,S,Gd")

def main():
    set_highlighted_excepthook()


    # Basic num arguments check
    assert len(sys.argv[1:]) >= 1, "Usage: python spatial_bound.py <sim_handle_1>"
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    data_folders = sys.argv[1:]
    for target in data_folders:
        label = target
        make_some_plots(target,molecular_path,label,dname_Figures,plot_derivative=PLOT_DERIVATIVE)

NUM_CELLS = 6 # TODO Automate
def make_some_plots(target,sim_output_parent_dir, label,figure_output_dir,plot_derivative = False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    for plot_mode in (PLOT_MODE,):#(0,1):     # 0: plot all charges, 1: plot element total charges

        fig, axs = plt.subplots(NUM_CELLS,1, sharex=True, facecolor='w')
        dashes = ["dashed","solid"]
        cmap = plt.get_cmap("Dark2")
        
        for s in range(NUM_CELLS):
            pl = Plotter(target,sim_output_parent_dir,spatial_index=s)
            pl.fig, pl.axs = fig,axs 
            pl.num_plotted=s  #hacky
            pl.plot_charges_bar("C",show_pulse_profile=False)
            
        plt.gcf().set_figwidth(FIGWIDTH)
        plt.gcf().set_figheight(PLOTHEIGHT*NUM_CELLS)
        #plt.gcf().tight_layout()
        #plt.tight_layout()
        qualifier = "charges_bar"
        plt.savefig(figure_output_dir + label + qualifier + figures_ext,bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option