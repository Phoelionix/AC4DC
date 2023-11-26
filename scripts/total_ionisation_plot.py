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

ylim=[None,3.2]
################
# stem = "SH_N"
# nums = range(1,12)
# ylim = [0,2] #[0,2.2]
# ylim = [0,4]#[0.25,2.25] #[0,2.2]
# stem = "SH_Xe"
# nums = range(0,6)
# ylim = [2,4]
# Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
stem_dict = {"SH_N":[1,11],
        #"SH_Zn":[1,11],
        "SH_Zr":[1,5],
        #"SH_Xe":[1,6],
    }
################
set_highlighted_excepthook()

def main():       
    molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"
    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    
    batches = {}
    # Get folders names in molecular_path
    all_outputs = os.listdir(molecular_path)
    valid_folder_names= True
    for key,val in stem_dict.items():
        data_folders = []
        for n in range(val[0],val[1]+1):
            data_folders.append(key+"-"+str(n)) # Folder name, excluding run tag "_"+R
        # Add run tag corresponding to latest run (highest "R" for "stem-n_R")
        for i, handle in enumerate(data_folders):
            matches = [match for match in all_outputs if match.startswith(handle)]
            run_nums = [int(R.split("_")[-1]) for R in matches if R.split("_")[-1].isdigit()]            
            # no folders - > not valid 
            if len(run_nums) == 0: 
                valid_folder_names = False
                print("\033[91mInput error\033[0m (a run corresponding to \033[91m"+handle+"\033[0m): not found in"+molecular_path)
                break
            data_folders[i] = handle + "_"+str(max(run_nums))
        batches[key] = data_folders
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    fig_title = "".join([stem.split("-")[0].split("_")[-1] for stem in batches.keys()])
    label = fig_title+'_total_ion_plot'
    plot(batches,molecular_path,label,dname_Figures,fig_title=fig_title)

def plot(batches,sim_output_parent_dir, label,figure_output_dir, fig_title = "", charge_conservation=False,bound_ionisation=False,free=False,free_slices=False,bound_ionisation_bar=False):
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

    fig, ax = plt.subplots(figsize=(FIGWIDTH,FIGHEIGHT))

    for stem, mol_names in batches.items():
        print("Plotting "+stem+" trace.")
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
        

        
        ax.scatter(energies,charges,label=stem.split("-")[0].split("_")[-1])
        if LABEL_TIMES:
            for i, txt in enumerate(times):
                ax.annotate(txt, (energies[i],charges[i]))
        ax.set_ylabel("Average carbon charge")
        ax.set_xlabel("Energy (keV)")
        
    ax.set_ylim(ylim)
    ax.legend()
    #ax.set_title(fig_title)
    plt.tight_layout()
    plt.savefig(figure_output_dir + label + figures_ext)
    plt.close()

if __name__ == "__main__":
    main()
