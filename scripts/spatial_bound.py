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
from core_functions import get_spatial_indices

#FIGWIDTH = 3.49751*2/3
#FIGHEIGHT = 3.49751/2
PLOTWIDTH = 3.49751*2/3
PLOTHEIGHT = 3.49751/2

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
        #make_some_plots(target,molecular_path,label,dname_Figures,plot_derivative=PLOT_DERIVATIVE)
        #overlaid_tot_charge(target,molecular_path,label,dname_Figures,intensity_averaged=True)
        overlaid_form_factor_disagreement(target,molecular_path,label,dname_Figures,intensity_averaged=False,q_list=np.linspace(0.1,4,10))

def make_some_plots(target,sim_output_parent_dir, label,figure_output_dir):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    spatial_indices = get_spatial_indices(sim_output_parent_dir,target)

    fig, axs = plt.subplots(len(spatial_indices),2, sharex=True, facecolor='w')
    dashes = ["dashed","solid"]
    cmap = plt.get_cmap("Dark2")
    
    ymax = None
    for s in spatial_indices:
        pl = Plotter(target,sim_output_parent_dir,spatial_index=s)
        pl.fig, pl.axs = fig,axs 
        pl.num_plotted=2*s  #hacky
        pl.plot_charges_bar("C",show_pulse_profile=False)
        #pl.plot_tot_charge(atoms=["C"], ylim=[0,1])
        ax = pl.plot_tot_charge(atoms=["C"], ylim=[0,ymax],intensity_averaged=True)
        if ymax is None:
            ymax = ax.get_ylim()[1]
        
    plt.gcf().set_figwidth(PLOTWIDTH*axs.shape[1])
    plt.gcf().set_figheight(PLOTHEIGHT*axs.shape[0])
    plt.gcf().tight_layout()
    #plt.tight_layout()
    qualifier = "charges_bar"
    plt.savefig(figure_output_dir + label + qualifier + figures_ext,bbox_inches='tight')
    plt.close()

def overlaid_tot_charge(target,sim_output_parent_dir, label,figure_output_dir,intensity_averaged):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    spatial_indices = get_spatial_indices(sim_output_parent_dir,target)

    fig, axs = plt.subplots(squeeze=False)
    dashes = ["dashed","solid"]
    cmap = plt.get_cmap("Dark2")
    
    ymax = None
    for s in spatial_indices:
        pl = Plotter(target,sim_output_parent_dir,spatial_index=s)
        pl.fig, pl.axs = fig,axs 
        pl.num_plotted=0
        pl.plot_tot_charge(atoms=["C"],intensity_averaged=intensity_averaged,label=str(s),plot_legend=False)
    #plt.legend(ncols=3)
    axs[0][0].legend(ncols=3)
        
    plt.gcf().set_figwidth(PLOTWIDTH*2)
    plt.gcf().set_figheight(PLOTHEIGHT*2)
    plt.gcf().tight_layout()
    #plt.tight_layout()
    qualifier = "charges"
    plt.savefig(figure_output_dir + label + qualifier + figures_ext,bbox_inches='tight')
    plt.close()

def overlaid_form_factor_disagreement(target,sim_output_parent_dir, label,figure_output_dir,intensity_averaged,q_list):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    ############
    # File/directory names
    #######  
    figures_ext = ".pdf"
    spatial_indices = get_spatial_indices(sim_output_parent_dir,target)

    # dashes = ["dashed","solid"]
    # cmap = plt.get_cmap("Dark2")
    
    for q in q_list:
        fig, axs = plt.subplots(squeeze=False)
        ymax = None
        for s in spatial_indices:
            pl = Plotter(target,sim_output_parent_dir,spatial_index=s)
            pl.fig, pl.axs = fig,axs 
            pl.num_plotted=0
            #pl.plot_form_factor_at_q(q,"C",intensity_averaged=intensity_averaged)
            pl.plot_form_factor_disagreement_at_q(q,"C",intensity_averaged=intensity_averaged,resolution=True,ylim=[None,0.1])
        #plt.legend(ncols=3)
        axs[0][0].legend(ncols=3)
            
        plt.gcf().set_figwidth(PLOTWIDTH*2)
        plt.gcf().set_figheight(PLOTHEIGHT*2)
        plt.gcf().tight_layout()
        #plt.tight_layout()
        qualifier = f"ff_disagreement-{q}"
        plt.savefig(figure_output_dir + label + qualifier + figures_ext,bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option