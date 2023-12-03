import matplotlib
matplotlib.use("pgf")
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    #thaumatin thing
    'axes.titlesize':10,     # fontsize of the axes title
    'axes.labelsize':10,    # fontsize of the x and y labels   
    'ytick.labelsize':10,
    'xtick.labelsize':10,
    'legend.fontsize':10,    
    'lines.linewidth':1,
    'lines.markersize':4,
})
import matplotlib.pyplot as plt
import numpy as np
from plotter_core import Plotter
from core_functions import get_sim_params
import sys, traceback
import os.path as path
import os
from QoL import set_highlighted_excepthook
import copy
import math
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
from scatter import XFEL,Crystal,stylin
from damage_landscape import multi_damage
from core_functions import get_pdb_path
import imaging_params
from scipy.interpolate import InterpolatedUnivariateSpline,SmoothBivariateSpline


FIGWIDTH_FACTOR = 1
FIGWIDTH = 3.49751 *FIGWIDTH_FACTOR
FIGHEIGHT = FIGWIDTH*3/4
FIT = True

matplotlib.rcParams.update({
    'figure.subplot.left':0.14/FIGWIDTH_FACTOR,
    'figure.subplot.right':1,
    'figure.subplot.bottom': 0.16/FIGWIDTH_FACTOR,
    'figure.subplot.top':1-0.07/(FIGWIDTH_FACTOR*FIGHEIGHT/FIGWIDTH*4/3),
})



LABEL_TIMES = False # Set True to graphically check times are all aligned.


MODE = 1 # | 0: mean carbon charge | 1: R factors |

ylim=[None,None]
################
# stem = "SH_N"
# nums = range(1,12)
# ylim = [0,2] #[0,2.2]
# ylim = [0,4]#[0.25,2.25] #[0,2.2]
# stem = "SH_Xe"
# nums = range(0,6)
# ylim = [2,4]
# Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
stem_dict = {
        "SH_N":[1,11],
        "SH_Zn":[1,11],
        #"SH_Zr":[1,10],
        "SH_Xe":[1,9],
    }
################
## Constants
# Dependent variable
AVERAGE_CHARGE = 0
R_FACTOR = 1
#
SCATTER_DIR = path.abspath(path.join(__file__ ,"../")) +"/scattering/"
PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"lys") 
MOLECULAR_PATH = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/" # directory of damage sim output folders

set_highlighted_excepthook()
def main():       
    plot_type = {AVERAGE_CHARGE:'avg_charge_plot',R_FACTOR:'R_dmg_plot'}        

    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    
    batches = {}
    # Get folders names in molecular_path
    all_outputs = os.listdir(MOLECULAR_PATH)
    valid_folder_names= True
    for key,val in stem_dict.items():
        data_folders = []
        for n in range(val[0],val[1]+1):
            data_folders.append(key+"-"+str(n)) # Folder name, excluding run tag "_"+R
        # Add run tag corresponding to latest run (highest "R" for "stem-n_R")
        for i, handle in enumerate(data_folders):
            matches = [match for match in all_outputs if match.startswith(handle+"_")]
            run_nums = [int(R.split("_")[-1]) for R in matches if R.split("_")[-1].isdigit()]            
            # no folders - > not valid 
            if len(run_nums) == 0: 
                valid_folder_names = False
                print("\033[91mInput error\033[0m (a run corresponding to \033[91m"+handle+"\033[0m): not found in"+MOLECULAR_PATH)
                break
            data_folders[i] = handle + "_"+str(max(run_nums))
        batches[key] = data_folders
    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder."
    fig_title = "".join([stem.split("-")[0].split("_")[-1] for stem in batches.keys()])
    label = fig_title+"_"+plot_type[MODE]
    plot(batches,label,dname_Figures,mode=MODE)

def plot(batches,label,figure_output_dir,mode = 0):
    ylabel = {AVERAGE_CHARGE:"Average carbon's charge",R_FACTOR:"$R_{dmg}$"}
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
        X = []
        Y = []
        print("Plotting "+stem+" trace.")
        energies = []
        times = []
        for mol_name in mol_names:
            energies.append(get_sim_params(mol_name)[2]/1000)
            pl = Plotter(mol_name)
            
            if mode is AVERAGE_CHARGE:
                intensity_averaged = False
                if intensity_averaged: 
                    Y.append(np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData))
                else:
                    step_index = -1  # Average charge at End-of-pulse 
                    Y.append(pl.get_total_charge(atoms="C")[step_index])
            elif mode is R_FACTOR:            
                Y.append(get_R(mol_name,MOLECULAR_PATH,imaging_params.goldilocks_dict)[0][0])

            times.append(pl.timeData[-1])
        X = energies
        print("Energies:",energies)
        if mode is AVERAGE_CHARGE:
            print("Charges:",Y)
        if mode is R_FACTOR:
            print("R factors:",Y)
        ax.scatter(X,Y,label=stem.split("-")[0].split("_")[-1])
        
        if FIT:
            #p = np.polynomial.Legendre.fit(X,Y,deg=6)
            #ax.plot(*p.linspace(100))
            sp = InterpolatedUnivariateSpline(X,Y)
            finer_x = np.linspace(np.min(X),np.max(X),endpoint=True)
            ax.plot(finer_x, sp(finer_x)) 

        if LABEL_TIMES:
            for i, txt in enumerate(times):
                ax.annotate(txt, (X[i],Y[i]))
        ax.set_ylabel(ylabel[mode])
        ax.set_xlabel("Energy (keV)")
    
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)
    #yticks = ax.get_yticks()
    #OOM = math.floor(math.log10(np.max(yticks)))
    #ax.set_yticks((yticks//(10**(OOM-1))*(10**(OOM-1))))
    #ax.set_yticklabels(([str(y/OOM) for y in yticks ]))

    ax.set_ylim(ylim)
    ax.legend(title="Dopant",fancybox=True)
    #ax.set_title(fig_title)
    #fig.tight_layout()
    fig.savefig(figure_output_dir + label + figures_ext)#bbox_inches="tight", pad_inches=0)

def grow_crystals(run_params,pdb_path,allowed_atoms,identical_deviations=True,plot_crystal=True,CNO_to_N=False,S_to_N=True):
    ''' 
    Grows a crystal. The first is damaged, corresponding to I_real, the second is undamaged, corresponding to I_ideal
    Ignores atoms not in allowed_atoms
    '''
    # The undamaged crystal form factor is constant (the initial state's) but still performs the same integration step with the pulse profile weighting.
    crystal_dmged = Crystal(PDB_STRUCTURE,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **run_params["crystal"])
    if identical_deviations:
        # we copy the other crystal so that it has the same deviations in coords
        crystal_undmged = copy.deepcopy(crystal_dmged)
        crystal_undmged.is_damaged = False
    else:
        crystal_undmged = Crystal(pdb_path,allowed_atoms,is_damaged=False,convert_excluded_elements_to_N=True, **run_params["crystal"])
    if plot_crystal:
        #crystal.plot_me(250000)
        crystal_undmged.plot_me(25000,template="plotly_dark")
    return crystal_dmged, crystal_undmged

def get_R(sim_handle,sim_handle_parent_folder,scattering_sim_parameters,allowed_atoms=["C","N","O"]):        
        run_params = copy.deepcopy(scattering_sim_parameters)
        print("Performing scattering simulation using data from damage simulation '" + sim_handle + "'.")
        exp1_qualifier = "_real"
        exp2_qualifier = "_ideal"          
        exp_name1 = sim_handle + "_" + exp1_qualifier
        exp_name2 = sim_handle + "_" + exp2_qualifier

        start_time,end_time,energy,fwhm,photon_count,param_dict,unit_dict = get_sim_params(sim_handle)

        # 
        print("Preparing XFEL...")
        experiment1 = XFEL(exp_name1,energy,**run_params["experiment"])
        experiment2 = XFEL(exp_name2,energy,**run_params["experiment"])

        print("Growing crystals at breakneck speeds...")     
        crystal_real, crystal_ideal = grow_crystals(run_params,PDB_STRUCTURE,allowed_atoms,plot_crystal=False)
    
        if run_params["laser"]["SPI"]:
            SPI_result1 = experiment1.spooky_laser(start_time,end_time,sim_handle,sim_handle_parent_folder,crystal_real, **run_params["laser"])
            SPI_result2 = experiment2.spooky_laser(start_time,end_time,sim_handle,sim_handle_parent_folder,crystal_ideal, **run_params["laser"])
            #TODO apply_background([SPI_result1,SPI_result2])
            R,cc,resolutions = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=True,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        '''
        else:
            exp1_orientations = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_real, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            if exp_name2 != None:
                # sync orientation with damage simulation
                experiment2.set_orientation_set(exp1_orientations)  
                run_params["laser"]["random_orientation"] = False 
                experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_ideal, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            R,cc = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=True,SPI=False,results_parent_dir = sctr_results_batch_dir)
        '''
        #pulse_params = [energy,fwhm,photon_count]
        return R,resolutions  # resolution index 0 corresponds to the max resolution

if __name__ == "__main__":
    main()

