# Comparing evolution of element average charges between multiple simulations.

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

OCCUPANCY = True # ONLY IMPLEMENTED FOR PLOT MODE 1 CURRENTLY
CHARGE_DIFFERENCE = False # Set True if want initially ionised species' traces to start from origin
PLOT_DERIVATIVE = False # Plot the rate of avg charge gain
#YLIM = [0,8]
PLOT_MODE = 1 # 1: plot element total charges # 2: plot orbital charges
#YLIM = [0,None]
YLIM = [3,6]
SCALE = 6
FIGWIDTH = SCALE*2/3
FIGHEIGHT = SCALE/2
#XLIM = [None,None]
XLIM = [None,None]
#YLIM=[0,6]
#YLIM=[20,24]

ATOM = "C"
#ATOMS = ("Cr_LDA",)

CUSTOM_LEGEND = None
#CUSTOM_LEGEND=("1s: CNO","2s: CNO","2p: CNO","1s: CNO,S,Gd","2s: CNO,S,Gd","2p: CNO,S,Gd")
#CUSTOM_LEGEND = ("Primary ionization only", "All ionization",)
def main():
    set_highlighted_excepthook()


    # Basic num arguments check
    assert len(sys.argv[1:]) > 0, "Usage: python compare_ion.py <sim_handle_1> <sim_handle_2> ..."
        
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
    label = data_folders[0]
    if len(data_folders) > 1:
         label+= "_"+data_folders[1]
    make_some_plots(data_folders,molecular_path,label,dname_Figures,plot_derivative=PLOT_DERIVATIVE)

def make_some_plots(mol_names,sim_output_parent_dir, label,figure_output_dir,plot_derivative = False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    for plot_mode in (PLOT_MODE,):
        fig, axs = plt.subplots(3, 3, sharey=True, facecolor='w')



        cmap = plt.get_cmap("Dark2")

        # Hacky way to get overlaid plots TODO
        pl = Plotter(mol_names[0],sim_output_parent_dir)
        pl.setup_axes(1)
        for m, mol_name in enumerate(mol_names):    
            pl.__init__(mol_name,sim_output_parent_dir)       
            pl.num_plotted = 0 # ƪ（˘へ˘ ƪ）
            if plot_mode == 1:
                colours = [cmap(m)]
                pl.plot_tot_charge(every=1,colours = colours,atoms = [ATOM],plot_legend=(m==0),xlim=XLIM,ylim=YLIM,charge_difference=CHARGE_DIFFERENCE,plot_derivative=PLOT_DERIVATIVE,occupancy=OCCUPANCY)

            if plot_mode == 2:
                ax = pl.plot_orbitals_charge(every=1,atom = ATOM,plot_legend=False,xlim=XLIM,ylim=YLIM,plot_derivative=PLOT_DERIVATIVE)       
            

            ax = pl.axs[0][0]
            if CUSTOM_LEGEND is None: 
                ax.legend(bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)  # Top right legend.
            else:
                handles,_ = ax.get_legend_handles_labels()
                handles = list(handles)
                assert len(CUSTOM_LEGEND)==len(handles)
                ax.legend(handles,CUSTOM_LEGEND,bbox_to_anchor=(1.02, 1),loc='upper left', ncol=1,handlelength=1)             
        plt.gcf().set_figwidth(FIGWIDTH)
        plt.gcf().set_figheight(FIGHEIGHT)
        #plt.gcf().tight_layout()
        #plt.tight_layout()

        if plot_mode == 1:
            qualifier = f"_{ATOM}-Comp"
        if plot_mode == 2:
            qualifier = f"_{ATOM}-OrbsComp"
        if PLOT_DERIVATIVE:
            qualifier+="-deriv"
        plt.savefig(figure_output_dir + label +qualifier + figures_ext,bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option