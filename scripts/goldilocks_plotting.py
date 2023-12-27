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
from scipy.interpolate import InterpolatedUnivariateSpline,SmoothBivariateSpline
import math
import json
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
from scatter import XFEL,Crystal,stylin
from damage_landscape import multi_damage
from core_functions import get_pdb_path
import imaging_params
from mergedeep import merge


FIGWIDTH_FACTOR = 1
FIGWIDTH = 3.49751 *FIGWIDTH_FACTOR
FIGHEIGHT = FIGWIDTH*3/4
FIT = False

matplotlib.rcParams.update({
    'figure.subplot.left':0.16/FIGWIDTH_FACTOR,
    'figure.subplot.right':1-0.02/FIGWIDTH_FACTOR,
    'figure.subplot.bottom': 0.16/FIGWIDTH_FACTOR,
    'figure.subplot.top':1-0.07/(FIGWIDTH_FACTOR*FIGHEIGHT/FIGWIDTH*4/3),
})


###
# INDEP_VARIABLE
PHOTON_ENERGY = 0; EDGE_SEPARATION = 1; ELECTRON_SOURCE_ENERGY = 2
# DEP_VARIABLE
AVERAGE_CHARGE = 0; R_FACTOR = 1
# BATCH
REAL_ELEMENTS = 0; ELECTRON_SOURCE = 1

####### User parameters #########
INDEP_VARIABLE = 0 # | 0: energy of photons |1: Energy separation from given edge |2: Artificial electron source energy
DEP_VARIABLE = 1 # | 0: mean carbon charge | 1: R factors (performs scattering simulations for each) |
BATCH = 0 #| 0: Real elements | 1: Electron source. 
SCATTERING_TARGET = 1 # Used if DEP_VARIABLE = 1.  | 0: unit lysozyme (light atoms) | 1: 3x3x3 lysozyme (light atoms no solvent)|

## Graphical
LEGEND = True
LABEL_TIMES = False # Set True to graphically check times are all aligned.

## Numerical

NORMALISING_STEM = None #"SH_N" # None  #Stem (str) to normalise traces by. (s.t. normalised trace becomes horizontal line). If None, does not normalise.
#NORMALISING_STEM = "SH_N"

ylim=[None,None]

if BATCH is REAL_ELEMENTS:
    assert(INDEP_VARIABLE != ELECTRON_SOURCE_ENERGY)
    # Batch stems and index range (inclusive). currently assuming form of key"-"+n+"_1", where n is a number in range of stem[key]
    # e.g. item "Carbon":[2,4] means to plot simulations corresponding to outputs Carbon-2, Carbon-3, Carbon-4, mol files. Last simulation output is used in case of duplicates (highest run number).
    stem_dict = {
            "SH_N":[2,10], #[700,704],
            #"SH_S":[2,10],  
            #"SH_Fe":[2,10],  
            #"SH_Zn":[2,10],
            #"SH_Se":[2,10],
            #"SH_Zr":[2,10],
            #-----L------
            #"SH_Ag":[2,8],
            #"SH_Xe":[2,9],
            "L_Gd":[1,9],
    }
    # Each item corresponds to a list of additional simulations to add on
    additional_points_dict = {

    }    
    
if BATCH is ELECTRON_SOURCE:
    assert(INDEP_VARIABLE == ELECTRON_SOURCE_ENERGY)
    ##
    #photon_energy = 10
    fluence = 1
    width = 15
    source_fraction = 3
    source_duration = 1
    ##
    # Manual mapping... TODO just pop outputs in folder or something.
    def batch_mapping(fluence,width,source_fraction,source_duration):
        params = fluence,width,source_fraction,source_duration
        if params == (1,15,1,1):
            return dict(ES_C = ([1,22],[0]))
        if params == (1,15,3,1):
            return  dict(ES_C = ([23,44],[0,902]))
        if params == (1,30,1,1):
            return dict(ES_C =([45,66],[]))
        raise KeyError("Invalid/missing batch params")
    stem_dict = {}
    additional_points_dict = {}
    for key, val in batch_mapping(fluence,width,source_fraction,source_duration).items():
        stem_dict[key] = val[0]
        additional_points_dict[key] = val[1]
        
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
    "ES_C":(0.3,"K"),
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
    "ES_C":0,
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
    "ES_C":2,
}

################
## Constants

SCATTERING_TARGET_DICT = {0: (imaging_params.goldilocks_dict_unit,"-unit"),1: (imaging_params.goldilocks_dict_3x3x3,"3x3x3"),99:(imaging_params.goldilocks_dict_unit,"-unit")}
#
SCATTER_DIR = path.abspath(path.join(__file__ ,"../")) +"/scattering/"
PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"lys") 
MOLECULAR_PATH = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/" # directory of damage sim output folders
GOLDI_RESULTS_PATH = path.abspath(path.join(__file__ ,"../goldilocks_saves")) + "/" 

xlim = [7,20]
if INDEP_VARIABLE is EDGE_SEPARATION:
    xlim[1] -= xlim[0]
    xlim[0] = 0 
if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
    #xlim = [None,None]
    xlim = [0,20]

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
        sim_tags = [*range(val[0],val[1]+1)]
        if key in additional_points_dict:
            for sim_tag in additional_points_dict[key]:
                sim_tags.append(sim_tag)
        for n in sim_tags:
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


    assert valid_folder_names, "One or more arguments (directory names) were not present in the output folder. See input error message(s)"
    fig_title = "".join([stem.split("-")[0].split("_")[-1] for stem in batches.keys()])
    label = fig_title+"_"+plot_type[DEP_VARIABLE]
    plot(batches,label,dname_Figures,mode=DEP_VARIABLE)

def plot(batches,label,figure_output_dir,mode = 0):
    ylabel = {AVERAGE_CHARGE:"Average carbon's charge",R_FACTOR:"$R_{dmg}$"}
    '''
    Arguments:
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
    if DEP_VARIABLE == 1:
        im_params,output_tag = SCATTERING_TARGET_DICT[SCATTERING_TARGET]
    if NORMALISING_STEM not in [None,False]:
        output_tag += "-normed"
    if INDEP_VARIABLE is EDGE_SEPARATION:
        output_tag +="-PhEl"  # X-ray photoelectron = XPhEl...
    if BATCH is ELECTRON_SOURCE:
        output_tag += str(fluence)+"-"+str(width)+"-"+str(source_fraction)+"-"+str(source_duration)

    if NORMALISING_STEM not in [None,False]:
        print("Getting norm")
        mol_names = batches[NORMALISING_STEM]
        norm = []
        for mol_name in mol_names:
            if INDEP_VARIABLE is PHOTON_ENERGY: 
                energy = get_sim_params(mol_name)[0]["energy"]
            if INDEP_VARIABLE is EDGE_SEPARATION:
                energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0]*1000
            if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
                energy = get_sim_params(mol_name)[0]["source_energy"]
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
            
            # Don't recalculate norm.
            indep_variable_key = str(INDEP_VARIABLE)
            dep_variable_key = str(DEP_VARIABLE)
            if DEP_VARIABLE == 1:
                dep_variable_key = str(DEP_VARIABLE)+"-"+str(SCATTERING_TARGET)
            if NORMALISING_STEM not in [None,False] and NORMALISING_STEM == stem:  #TODO tidy up by moving outside loop
                for elem in norm:
                    X.append(elem[0])
                    Y.append(elem[1])
                break
            else:
                saved_x = saved_y = saved_end_time = None
                dat= get_saved_data(mol_name) 
                if dat is not None:
                        saved_x,                              saved_y,                            saved_end_time = (
                        dat.get("x").get(indep_variable_key),dat.get("y").get(dep_variable_key),dat.get("end_time")
                    )
                ## Get independent variable
                if INDEP_VARIABLE is PHOTON_ENERGY: 
                    energy = get_sim_params(mol_name)[0]["energy"]
                if INDEP_VARIABLE is EDGE_SEPARATION:
                    energy = get_sim_params(mol_name)[0]["energy"] - edge_dict[stem][0]*1000
                if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
                    s = get_sim_params(mol_name)[0]
                    energy = s["source_energy"]
                    if energy is None:
                        # No source
                        if (s["fluence"],s["width"]) != (fluence,width):
                            print(s)
                            print((s["fluence"],s["width"]),"!=",(fluence,width))
                            raise ValueError("Params do not match expected values for batch")                        
                        # Plot horizontal dashed line corresponding to damage with no artificial source.
                        ax.axhline(y=Y[-1],ls=(5,(10,3)),color=cmap(0),label="No injected")
                        Y.pop(-1)
                        continue
                    if (s["fluence"],s["width"],s["source_fraction"],s["source_duration"]) != (fluence,width,source_fraction,source_duration):
                        print(s)
                        print((s["fluence"],s["source_fraction"],s["source_duration"]),"!=",(fluence,width,source_fraction,source_duration))
                        raise ValueError("Params do not match expected values for batch")                    
                X.append(energy/1000)    
                if saved_x is not None:
                    assert saved_x == X[-1]                            
                ## Measure damage 
                if saved_y is None:              
                    pl = Plotter(mol_name)
                    if mode is AVERAGE_CHARGE:
                        intensity_averaged = False
                        if intensity_averaged: 
                            Y.append(np.average(pl.get_total_charge(atoms="C")*pl.intensityData)/np.average(pl.intensityData))
                        else:
                            step_index = -1  # Average charge at End-of-pulse 
                            Y.append(pl.get_total_charge(atoms="C")[step_index])
                    elif mode is R_FACTOR:            
                        Y.append(get_R(mol_name,MOLECULAR_PATH,im_params)[0][0])
                else: 
                    Y.append(saved_y)
                if saved_end_time is None:
                    if saved_y is not None:
                        pl = Plotter(mol_name)
                    times.append(pl.timeData[-1]) 
                else:
                    times.append(saved_end_time)                   

                save_dict = {}
                if saved_x is None:
                    save_dict["x"] = {indep_variable_key:X[-1]}
                if saved_y is None:    
                    save_dict["y"] = {dep_variable_key:Y[-1]}
                if saved_end_time is None:
                    save_dict["end_time"] = times[-1]
                if save_dict is not {}:
                    save_data(mol_name,save_dict)
        print("Energies:",X)
        if mode is AVERAGE_CHARGE:
            print("Charges:",Y)
        if mode is R_FACTOR:
            print("R factors:",Y)
        if NORMALISING_STEM not in [None,False]:
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

        if BATCH is REAL_ELEMENTS:
            _label = dopant 
            if ground_charge_dict[stem]>1: 
                _label+="$^{"+str(ground_charge_dict[stem])+"+}$"
            elif ground_charge_dict[stem]==1:
                _label+="$^{+}$"
        if BATCH is ELECTRON_SOURCE:
            _label = "Injected"
        ax.scatter(X,Y,label=_label,color = cmap(c))
        
        k =2 # Spline order
        if FIT and len(X)>=k+1:
            ordered_dat = sorted(zip(X,Y))
            # Split dataset with ionisation edge of dopant
            split_idx = 0
            if INDEP_VARIABLE is PHOTON_ENERGY:
                split_idx = np.searchsorted(np.array([x[0] for x in ordered_dat]),edge_dict[stem][0])
            if INDEP_VARIABLE is EDGE_SEPARATION:
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
        if NORMALISING_STEM not in [None,False]:
            _ylab += " relative difference" # TODO better clarity
        ax.set_ylabel(_ylab)
        if INDEP_VARIABLE is PHOTON_ENERGY:
            ax.set_xlabel("Photon energy (keV)")
        if INDEP_VARIABLE is EDGE_SEPARATION:
            ax.set_xlabel("Dopant edge separation (keV)")  # TODO define better
        if INDEP_VARIABLE is ELECTRON_SOURCE_ENERGY:
            ax.set_xlabel("Source electron energy (keV)")  # TODO define better
        
    if xlim[0] is None:
        xlim[0] = np.min(X)
    if xlim[1] is None:
        xlim[1] = np.max(X)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)           

    if INDEP_VARIABLE is PHOTON_ENERGY: 
        for stem, mol_names in batches.items():
            c = col_dict[stem]
            dopant = stem.split("-")[0].split("_")[-1]
            if ax.get_xlim()[0] < edge_dict[stem][0] < ax.get_xlim()[1]:
                ax.axvline(x=edge_dict[stem][0],ls=(5,(10,3)),color=cmap(c))
                ax.text(edge_dict[stem][0]-0.1, 0.96*ax.get_ylim()[1] +0.04*ax.get_ylim()[0],dopant+"--$"+edge_dict[stem][1]+"$",verticalalignment='top',horizontalalignment='right',rotation=-90,color=cmap(c))
    if INDEP_VARIABLE is EDGE_SEPARATION:
        ax.axvline(x=0,ls=(5,(10,3)),color="black")

    
    ax.ticklabel_format(style='sci', axis='y', scilimits=(-1,1),)
    #yticks = ax.get_yticks()
    #OOM = math.floor(math.log10(np.max(yticks)))
    #ax.set_yticks((yticks//(10**(OOM-1))*(10**(OOM-1))))
    #ax.set_yticklabels(([str(y/OOM) for y in yticks ]))

    if LEGEND:
        if BATCH is REAL_ELEMENTS:
            ax.legend(title="Dopant",fancybox=True,ncol=2,loc='upper center',bbox_to_anchor=(0.75, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)
        if BATCH is ELECTRON_SOURCE:
            ax.legend(fancybox=True,ncol=1,loc='upper center',bbox_to_anchor=(0.8, 1.02),handletextpad=0.01,columnspacing=0.2,borderpad = 0.18)

    #ax.set_title(fig_title)
    #fig.tight_layout()
    fig.savefig(figure_output_dir + label + output_tag + figures_ext)#bbox_inches="tight", pad_inches=0)

def grow_crystals(run_params,pdb_path,allowed_atoms,identical_deviations=True,plot_crystal=True,CNO_to_N=False,S_to_N=True):
    ''' 
    Grows a crystal. The first is damaged, corresponding to I_real, the second is undamaged, corresponding to I_ideal
    Ignores atoms not in allowed_atoms
    '''
    # The undamaged crystal form factor is constant (the initial state's) but still performs the same integration step with the pulse profile weighting.
    crystal_dmged = Crystal(pdb_path,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **run_params["crystal"])
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

        param_dict,_,_ = get_sim_params(sim_handle)
        start_time = param_dict["start_t"]
        end_time = param_dict["end_t"]
        energy = param_dict["energy"]
        # 
        print("Preparing XFEL...")
        experiment1 = XFEL(exp_name1,energy,**run_params["experiment"])
        experiment2 = XFEL(exp_name2,energy,**run_params["experiment"])

        print("Growing crystals at breakneck speeds...")     
        structure = PDB_STRUCTURE
        if SCATTERING_TARGET == 99:
            # Debug, use tetrapeptide (small structure, fast)
            structure = get_pdb_path(SCATTER_DIR,"tetra")
            
        crystal_real, crystal_ideal = grow_crystals(run_params,structure,allowed_atoms,plot_crystal=False)
    
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

def get_saved_data(handle):
    # Look for save of data
    file_path = GOLDI_RESULTS_PATH+handle+".json"
    if os.path.exists(file_path):
        with open(file_path,"r") as f:
            return json.load(f)
    return None
    
    
def save_data(handle,output_data_dict):
    # Delete saves for smaller older simulations with this mol file.
    all_outputs = os.listdir(MOLECULAR_PATH)
    for match in all_outputs:
        if match.startswith(handle+"_") and match != handle:
            os.remove(GOLDI_RESULTS_PATH+match+".json")

    os.makedirs(GOLDI_RESULTS_PATH,exist_ok=True)
    file_path = GOLDI_RESULTS_PATH+handle+".json"
    # Get any data from existing json file 
    data = output_data_dict
    if os.path.exists(file_path):
        with open(GOLDI_RESULTS_PATH+handle+".json","r") as f:
            saved_data_dict = json.load(f)
            data = merge(data,saved_data_dict)

    # Save data
    with open(GOLDI_RESULTS_PATH+handle+".json", 'w') as f: 
        json.dump(data, f)    

if __name__ == "__main__":
    main()

