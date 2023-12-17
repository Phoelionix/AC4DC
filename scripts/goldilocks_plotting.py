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
    'figure.subplot.right':1-0.02/FIGWIDTH_FACTOR,
    'figure.subplot.bottom': 0.16/FIGWIDTH_FACTOR,
    'figure.subplot.top':1-0.07/(FIGWIDTH_FACTOR*FIGHEIGHT/FIGWIDTH*4/3),
})


####### User parameters #########
## Graphical
LEGEND = True
LABEL_TIMES = False # Set True to graphically check times are all aligned.

## Numerical
MODE = 1 # | 0: mean carbon charge | 1: R factors (performs scattering simulations for each) |
SCATTERING_TARGET = 0 # Used if MODE = 1.  | 0: unit lysozyme (light atoms) | 1: 3x3x3 lysozyme (light atoms no solvent)|
NORMALISING_STEM = None#"SH_N" # None  # None, or the stem (str) to normalise traces by. (s.t. normalised trace becomes horizontal line)
INDEP_VARIABLE = 2 # | 0: energy of photons |1: Energy separation from given edge |2: Artificial electron source energy

ylim=[None,None]
xlim = [None,None]

#ylim = [,]
xlim = [0, None]
#xlim = [7, 18]#[7,18] # [5,18] [7,18]

# Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
REAL_ELEMENTS = 0; ELECTRON_SOURCE = 1
#batch = REAL_ELEMENTS
batch = ELECTRON_SOURCE
if batch is REAL_ELEMENTS:
    # e.g. item "Carbon":[2,4] means to plot simulations corresponding to outputs Carbon-2, Carbon-3, Carbon-4, mol files. Last simulation output is used in case of duplicates (highest run number).
    stem_dict = {
            "SH_N":[2,11],
            "SH_S":[2,11],  
            "SH_Fe":[2,9],  
            "SH_Zn":[2,11],
            "SH_Se":[2,11],
            "SH_Zr":[2,10],
            #-----L------
            "SH_Ag":[2,8],
            "SH_Xe":[2,9],
            #"L_Gd":[1,1],
    }
    # Each item corresponds to a list of additional simulations to add on
    additional_points_dict = {

    }    
if batch is ELECTRON_SOURCE:
    stem_dict = {
        #"ES":[32,35],
        "ES_L":[47,52],
    }
    additional_points_dict = {
        "ES_L":[0,],
    }
# Ionisation edges of dopants
edge_dict = {
    "SH_N": (0.3,"K"),
    "SH_S": (1,"K"),
    "SH_Fe":[7.1,"K"],  
    "SH_Zn":(9.7,"K"),
    "SH_Se":(12.7,"K"),
    "SH_Zr":(17.98,"K"),  # Slight offset for graphical purposes.
    "SH_Ag":(3.8,"L_{1}"),
    "SH_Xe":(5.5,"L_{1}"),
    "L_Gd":(7.4,"L_{1}"),
    
    "ES":(0.3,"K"),
    "ES_L":(0.3,"K"),
}
ground_charge_dict = {
    "SH_N": 0,
    "SH_S": 0,
    "SH_Fe":2,  
    "SH_Zn": 2,
    "SH_Se":0,
    "SH_Zr":2,
    "SH_Ag":1,
    "SH_Xe":8,
    "L_Gd":11,

    "ES":0,
    "ES_L":0,
}
col_dict = {
    "SH_N": 0,
    "SH_S": 1,
    "SH_Fe":7,  
    "SH_Zn": 3,
    "SH_Se":4,
    "SH_Zr":5,
    "SH_Ag":8,
    "SH_Xe":6,
    "L_Gd":2,    

    "ES":0,
    "ES_L":1,
}

################
## Constants
# Dependent variable
# MODE
AVERAGE_CHARGE = 0
R_FACTOR = 1
# INDEP_VARIABLE
PHOTON_ENERGY = 0
PHOTOELECTRON_ENERGY = 1  # Not quite an accurate label.
ELECTRON_SOURCE_ENERGY = 2

SCATTERING_TARGET_DICT = {0: (imaging_params.goldilocks_dict_unit,"-unit"),1: (imaging_params.goldilocks_dict_3x3x3,"3x3x3")}
#
SCATTER_DIR = path.abspath(path.join(__file__ ,"../")) +"/scattering/"
PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"lys") 
MOLECULAR_PATH = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/" # directory of damage sim output folders

set_highlighted_excepthook()
def main():       


    assert all([k in edge_dict.keys() and k in ground_charge_dict.keys() for k in stem_dict.keys()])
    # Plot tag
    plot_type = {AVERAGE_CHARGE:'avg_charge_plot',R_FACTOR:'R_dmg'}      

    dname_Figures = "../../output/_Graphs/plots/"
    dname_Figures = path.abspath(path.join(__file__ ,dname_Figures)) + "/"
    
    batches = {}
    # Get folders names in molecular_path
    all_outputs = os.listdir(MOLECULAR_PATH)
    valid_folder_names= True
    for key,val in stem_dict.items():
        data_folders = []
        # Get folder name, excluding run tag
        for n in range(val[0],val[1]+1):
            data_folders.append(key+"-"+str(n)) 
        # Add run tag "_"+n corresponding to latest run (n = highest "R" for "stem-n_R")
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
    for key, sims in additional_points_dict.items():
        for sim_tag in sims:
            handle = key+"-"+str(sim_tag)
            matches = [match for match in all_outputs if match.startswith(handle+"_")]
            run_nums = [int(R.split("_")[-1]) for R in matches if R.split("_")[-1].isdigit()] 
            if len(run_nums) == 0: 
                valid_folder_names = False
                print("\033[91mInput error\033[0m (a run corresponding to \033[91m"+handle+"\033[0m): not found in"+MOLECULAR_PATH)
                break
            data_folder = handle + "_"+str(max(run_nums))
            assert(data_folder not in data_folders)
            data_folders.append(data_folder)             


    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder. See input error message(s)"
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
    cmap = plt.get_cmap('tab10')

    output_tag = ""
    if MODE == 1:
        im_params,output_tag = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    if NORMALISING_STEM is not None:
        output_tag += "-normed"
    if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
        output_tag +="-PhEl"  # X-ray photoelectron = XPhEl...

    if NORMALISING_STEM is not None:
        print("Getting norm")
        mol_names = batches[NORMALISING_STEM]
        norm = []
        for mol_name in mol_names:
            if INDEP_VARIABLE is PHOTON_ENERGY: 
                energy = get_sim_params(mol_name)[2]/1000
            if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
                energy = get_sim_params(mol_name)[2]/1000 - edge_dict[stem][0]
            if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
                energy = get_sim_params(mol_name,get_source_energy=True)[2]
                if energy is None:
                    energy = 0
                energy/=1000
            pl = Plotter(mol_name)
            if mode is AVERAGE_CHARGE:
                intensity_averaged = False
                if intensity_averaged: 
                    norm.append(np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData))
                else:
                    step_index = -1  # Average charge at End-of-pulse 
                    norm.append((energy,pl.get_total_charge(atoms="C")[step_index]))
            elif mode is R_FACTOR:            
                norm.append((energy,get_R(mol_name,MOLECULAR_PATH,im_params)[0][0]))


    for stem, mol_names in batches.items():
        print("Plotting "+stem+" trace.")
        c = col_dict[stem]
        dopant = stem.split("-")[0].split("_")[-1]
        X = []
        Y = []
        times = []  # last times of sims, for checking everything lines up.
        for mol_name in mol_names:
            pl = Plotter(mol_name)
            times.append(pl.timeData[-1])

            if NORMALISING_STEM is not None and NORMALISING_STEM == stem:  #TODO tidy up by moving outside loop
                for elem in norm:
                    X.append(elem[0])
                    Y.append(elem[1])
                break
            else:
                if INDEP_VARIABLE is PHOTON_ENERGY: 
                    X.append(get_sim_params(mol_name)[2]/1000)
                if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
                    X.append(get_sim_params(mol_name)[2]/1000 - edge_dict[stem][0])
                if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
                    energy = get_sim_params(mol_name,get_source_energy=True)[2]
                    if energy is None:
                        energy = 0
                    energy/=1000
                    X.append(energy)

                if mode is AVERAGE_CHARGE:
                    intensity_averaged = False
                    if intensity_averaged: 
                        Y.append(np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData))
                    else:
                        step_index = -1  # Average charge at End-of-pulse 
                        Y.append(pl.get_total_charge(atoms="C")[step_index])
                elif mode is R_FACTOR:            
                    Y.append(get_R(mol_name,MOLECULAR_PATH,im_params)[0][0])
        print("Energies:",X)
        if mode is AVERAGE_CHARGE:
            print("Charges:",Y)
        if mode is R_FACTOR:
            print("R factors:",Y)
        if NORMALISING_STEM is not None:
            for i, energy in enumerate(X): 
                for elem in norm:
                    if elem[0] == energy:
                        Y[i]/=elem[1]
                        break
                assert(False,"Missing norm for energy") #TODO: If not found, remove corresponding value. (go backwards in reverse )
                #     del(X[i])
                #     del(Y[i])
                #     del(energies[i])
                #     del(times[i])

        _label = dopant 
        if ground_charge_dict[stem]>1: 
            _label+="$^{"+str(ground_charge_dict[stem])+"+}$"
        elif ground_charge_dict[stem]==1:
            _label+="$^{+}$"
        ax.scatter(X,Y,label=_label,color = cmap(c))
        
        k =2 # Spline order
        if FIT and len(X)>=k+1:
            ordered_dat = sorted(zip(X,Y))
            # Split dataset with ionisation edge of dopant
            split_idx = 0
            if INDEP_VARIABLE is PHOTON_ENERGY:
                split_idx = np.searchsorted(np.array([x[0] for x in ordered_dat]),edge_dict[stem][0])
            if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
                split_idx = np.searchsorted(np.array([x[0] for x in ordered_dat]),-0.000000001)
            splines = []
            finer_x = []
            if split_idx in (0, len(X)):
                # No split in fit
                splines.append(InterpolatedUnivariateSpline(*zip(*ordered_dat),k=k))
                finer_x.append(np.linspace(np.min(X),np.max(X),endpoint=True))
            else:
                # Split in fit                  \\START           \\END
                for start_idx,end_idx  in zip([0,split_idx],[split_idx-1,len(X)-1]):
                    if end_idx+1 - start_idx>=k+1:
                        splines.append(InterpolatedUnivariateSpline(*zip(*ordered_dat[start_idx:end_idx+1]),k=k))
                        finer_x.append(np.linspace(ordered_dat[start_idx][0],ordered_dat[end_idx][0],endpoint=True))
            for sp, lil_x in zip(splines,finer_x):
                ax.plot(lil_x, sp(lil_x),color = cmap(c)) 

        if LABEL_TIMES:
            for i, txt in enumerate(times):
                ax.annotate(txt, (X[i],Y[i]))
        _ylab = ylabel[mode]
        if NORMALISING_STEM is not None:
            _ylab += " relative difference" # TODO better clarity
        ax.set_ylabel(_ylab)
        if INDEP_VARIABLE is PHOTON_ENERGY:
            ax.set_xlabel("Photon energy (keV)")
        if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
            ax.set_xlabel("Dopant `photoelectron energy' (keV)")  # TODO define better
        if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
            ax.set_xlabel("Source electron energy (keV)")  # TODO define better
        
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)           
    if INDEP_VARIABLE is PHOTON_ENERGY: 
        for stem, mol_names in batches.items():
            c = col_dict[stem]
            dopant = stem.split("-")[0].split("_")[-1]
            if ax.get_xlim()[0] < edge_dict[stem][0] < ax.get_xlim()[1]:
                ax.axvline(x=edge_dict[stem][0],ls=(5,(10,3)),color=cmap(c))
                ax.text(edge_dict[stem][0]-0.1, 0.96*ax.get_ylim()[1] +0.04*ax.get_ylim()[0],dopant+"--$"+edge_dict[stem][1]+"$",verticalalignment='top',horizontalalignment='right',rotation=-90,color=cmap(c))
    if INDEP_VARIABLE is PHOTOELECTRON_ENERGY:
        ax.axvline(x=0,ls=(5,(10,3)),color="black")

    
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)
    #yticks = ax.get_yticks()
    #OOM = math.floor(math.log10(np.max(yticks)))
    #ax.set_yticks((yticks//(10**(OOM-1))*(10**(OOM-1))))
    #ax.set_yticklabels(([str(y/OOM) for y in yticks ]))

    if LEGEND:
        ax.legend(title="Dopant",fancybox=True,ncol=2,loc='upper center',bbox_to_anchor=(0.75, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)
    #ax.set_title(fig_title)
    #fig.tight_layout()
    fig.savefig(figure_output_dir + label + output_tag + figures_ext)#bbox_inches="tight", pad_inches=0)

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

