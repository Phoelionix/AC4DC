import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    #thaumatin thing
    'axes.titlesize':12,     # fontsize of the axes title
    'axes.labelsize':12,    # fontsize of the x and y labels   
    'ytick.labelsize':12,
    'xtick.labelsize':12,
    'legend.fontsize':12,    
    'lines.linewidth':1,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
from core_functions import get_sim_params
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook

FIGWIDTH = 3.49697
FIGHEIGHT = FIGWIDTH*3/4
LABEL_TIMES = False # Set True to graphically check times are all aligned.

ylim=[None,None]
################
# stem = "SH_N"
# nums = range(1,12)
# ylim = [0,2] #[0,2.2]
stem = "SH_Zn"
nums = range(1,12)
ylim = [0,2]#[0.25,2.25] #[0,2.2]
# stem = "SH_Xe"
# nums = range(0,5)
# ylim = [2,4]
################
set_highlighted_excepthook()

def main():
    data_folders = []
    for val in nums:
        data_folders.append(stem+"-"+str(val)+"_1") 
        
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    valid_folder_names= True
    for i, data_folder in enumerate(data_folders):
        if not path.isdir(molecular_path+data_folder):
            valid_folder_names = False
            print("\033[91mInput error\033[0m (argument \033[91m"+str(i)+ "\033[0m): folder name not found.")
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    fig_title =data_folders[0].split("-")[0].split("_")[-1] 
    label = fig_title+'_total_ion_plot'
    plot(data_folders,molecular_path,label,dname_Figures,fig_title=fig_title)

def plot(mol_names,sim_output_parent_dir, label,figure_output_dir, fig_title = "", charge_conservation=False,bound_ionisation=False,free=False,free_slices=False,bound_ionisation_bar=False):
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

    step_index = -1
    charges = []
    energies = []
    times = []
    for mol_name in mol_names:
        energies.append(get_sim_params(mol_name)[2]/1000)
        pl = Plotter(mol_name)
        charges.append(pl.get_total_charge(atoms="C")[step_index])
        times.append(pl.timeData[step_index])
    print("Energies:",energies)
    print("Charges:",charges)
    fig, ax = plt.subplots(figsize=(FIGWIDTH,FIGHEIGHT))

    
    ax.scatter(energies,charges)
    if LABEL_TIMES:
        for i, txt in enumerate(times):
            ax.annotate(txt, (energies[i],charges[i]))
    ax.set_ylabel("Average carbon charge")
    ax.set_xlabel("Energy (keV)")
    ax.set_title(fig_title)
    
    ax.set_ylim(ylim)

    plt.tight_layout()
    plt.savefig(figure_output_dir + label + figures_ext)
    plt.close()

if __name__ == "__main__":
    main()
