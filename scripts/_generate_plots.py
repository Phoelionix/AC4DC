import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
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

ELECTRON_DENSITY = False # energy density if False

def main():
    set_highlighted_excepthook()


    # Basic num arguments check
    if  len(sys.argv) < 2:
        print("Usage: python3 generate_plots.py Carbon_1")
        print({"Pass extra arguments to plot for each"})
        exit()
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(sys.argv[1:]):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    for data_folder in sys.argv[1:]:
        label = data_folder +'_Plt'
        make_some_plots(data_folder,molecular_path,label,dname_Figures,True,True,True,True)

def make_some_plots(mol_name,sim_output_parent_dir, label,figure_output_dir, charge_conservation=False,bound_ionisation=False,free=False,free_slices=False):
    '''
    Arguments:
    mol_name: The name of the folder containing the simulation's data (the csv files). (By default this is the stem of the mol file.)
    sim_data_parent_dir: absolute path to the folder containing the folders specified by target_handles.
    '''    
    
    ############
    # File/directory names
    #######  
    figures_ext = "" #.png
    fname_charge_conservation = "charge_conservation"
    fname_free = "free"
    fname_HR_style = "HR_style"
    fname_bound_dynamics = "bound_dynamics"

    pl = Plotter(mol_name,sim_output_parent_dir,use_electron_density = ELECTRON_DENSITY)
    num_atoms = len(pl.statedict)
    num_subplots = charge_conservation + bound_ionisation*num_atoms + free + free_slices
    pl.setup_axes(num_subplots)

    if charge_conservation: 
        pl.plot_tot_charge(every=10)
        #plt.savefig(figure_output_dir + label + fname_charge_conservation + figures_ext)

    if bound_ionisation:
        pl.plot_all_charges()
        pl.fig.subplots_adjust(left=0.2,bottom=0.18,top=0.95)
        #plt.savefig(figure_output_dir + label + fname_bound_dynamics + figures_ext)
        
    if free:
        pl.plot_free(log=True, min=1e-7, every=5)
        #plt.savefig(figure_output_dir + label + fname_free + figures_ext)

    if free_slices:
        pl.initialise_step_slices_ax()
        from plotter_core import fit_maxwell, maxwell

        cmap = plt.get_cmap("tab10")

        plt.rcParams["font.size"] = 10 # 8

        ### defaults
        #Cutoff energies for fitting MB curves. TODO get from AC4DC
        thermal_cutoff_energies = [2000]*10     
        xmin,xmax = 1,1e4
        
        ###

        #TODO get slices from AC4DC or user input           
        # Hau-Riege (Sanders results)
        thermal_cutoff_energies = [200, 500, 500, 1000]
        slices = [-7.5,-5,-2.5,0]         
        # # -7.5 fs Hau-Riege
        # thermal_cutoff_energies = [200]
        # slices = [-7.5] 

        # # Royle Sect. B
        # thermal_cutoff_energies = [1500,1500,1500,1500]
        # slices = [-90, -50, 0]  
        # # Royle Sect. C
        # thermal_cutoff_energies = [500,500,600,1000]
        # slices = [-15, 0, 15, 30]          

        colrs = [cmap(i) for i in range(len(slices))]

        plot_legend = True        
        plot_fits = True # Whether to fit MB curves to distribution below thermal cutoff energies.
        ####### 
        # Here we can plot MB curves e.g. for fit comparison
        ####
        plot_custom_fits = False
        # # example
        # custom_T = [30.8,75.5,131.8,207.5]
        # custom_n = [0.06*3/2,0.06*3/2,0.06*3/2,0.06*3/2]
        # #Sanders/H-R for -7.5 fs  
        # custom_T = [44.1,31]  
        # custom_n = [0.06*3/2,0.06*3/2] # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = ['r','b']
        #H-R
        # custom_T = [31,70,125,195]   # [44.1,84.9,135.6,205.8]
        # custom_n = [0.06*3/2]*len(custom_T) # - (not anything meaningful, just density of MB approximated as 50% of total). 
        # custom_colrs = colrs   
        # xmin,xmax = 1,1500 
        # plot_legend = False    
        # #Royle B
        # custom_T = [8,33,105]  
        # custom_n = [0.08,0.08,0.08]     
        # xmin, xmax = 1,2000
        # colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        # custom_colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]  
        # #Royle C
        # custom_T = [100]  
        # custom_n = [0.155]     
        # colrs = [cmap(0),cmap(2),cmap(1),cmap(3)]          
        # custom_colrs = ['black']
        # plot_fits = True  # our fits don't work as well for some reason.
        #######
        
        lw = 1.5
        for (t, e, col ) in zip(slices, thermal_cutoff_energies, colrs):
                lines = pl.plot_step(t, normed=True, color = col, lw=lw)
                if plot_fits:
                    T = pl.plot_fit(t, e, normed=True, color=col, lw=lw,alpha=0.7)
        if plot_custom_fits:
            for (T, n, col ) in zip(custom_T, custom_n, custom_colrs):
                pl.ax_steps.plot([0],[0],alpha=0,label=" ")
                pl.plot_maxwell(T,n,color = col, lw=lw,alpha=0.7)


        pl.ax_steps.set_ylim([1e-4, 1])
        if ELECTRON_DENSITY:
            pl.ax_steps.set_ylim([2e-7, 1e-2]) #royle sect. B
            pl.ax_steps.set_ylim([1e-8, 1e-3]) #royle sect. C
            pl.ax_steps.set_xscale("linear")

        #TODO get from AC4DC
        pl.ax_steps.set_xlim([xmin,xmax]) #Hau-Riege

        #pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.95)
        pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
        pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False

        if plot_legend:
            handles, labels = pl.ax_steps.get_legend_handles_labels()
            # e.g. for 8 labels (4 time steps), order = [0,2,4,6,1,3,5,7]  
            order = list(range(0,len(labels) - 1,2)) + list(range(1,len(labels),2))
            if len(labels)%2 != 0:
                order.append(len(order))  # shouldnt happen though.

        
            pl.ax_steps.legend([handles[idx] for idx in order],[labels[idx] for idx in order], loc='upper center',ncol=2)

        name = label.replace('_',' ')
        pl.ax_steps.set_title(name + " - Free-electron distribution")

        #plt.savefig(figure_output_dir + label + fname_HR_style + figures_ext)
    plt.tight_layout()
    plt.savefig(figure_output_dir + label + figures_ext)
    plt.close()

if __name__ == "__main__":
    main()

#TODO: Change size to match my screen by default, add --y option