#%%
from scatter import *
from multiprocessing import Pool
from sample_scalepack import all_reflections_to_scalepack
import inspect

#TODO auto generate ideal, undamaged, damaged. not undamaged and damaged. (undamaged is called ideal)

NUM_PARALLEL=1
SEEDED = False
RANDOM_WATER = False; NUM_RANDOM_WATER = 0 # 702
DEBUG_WATER = False
QUICK_TEST = False
SKIP_UNDAMAGED = False

TARGET_HANDLE_DICT = dict(
    #lys_salt = "lys_salt_solvated_fast_H_5",
    #lys_no_salt =  "lys_solvated_fast_H_6"
    lys_salt = "lys_salt_solvated_H_40_2",
    lys_no_salt = "lys_solvated_H_40_3",
    lys_high_damage= "lys_galli_HF_23",
)
target_options = ["lys_no_salt"]
TAG = "one"


def main(par_idx):

    cycles_per_bragg_set = 1 # leave at 1, should be fine
    num_bragg_sets = 1 # 25 # increase this for better stochastic
    num_unique_supercells = 1 # 20

    fig_width = 3.49751 # 20
    fig_height = fig_width*3/4 # 20
    ### Simulate
    #target_options = ["lys_salt","lys_no_salt","neutze","hen","tetra","glycine","fcc"]
    #target_options = ["lys_salt","lys_no_salt"]
    
    target = target_options[0]
    if par_idx >= NUM_PARALLEL/2 and len(target_options)==2:
        target = target_options[1]
    if NUM_PARALLEL==1:
        par_idx_for_target=1
    else:
        par_idx_for_target = par_idx% int((NUM_PARALLEL/2))+1
    #============------------User params---------==========#
    assert target in target_options
    #target = "lys_no_salt"#"glycine"  #target_options[2]
    best_resolution = 1.3 # 1.58 (abdullah) # 2   # resolution (determining max q)
    worst_resolution = None#30 # 'resolution' corresponding to min q

    #---------------------------------#
    water_index = None # None TODO automate
    if target in TARGET_HANDLE_DICT:


        include_symmetries=True
        if QUICK_TEST:
            num_bragg_sets = 1
            num_unique_supercells = 1
            include_symmetries=False
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_1.4.hkl"
        unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_2.0.hkl"

        target_handle = TARGET_HANDLE_DICT[target]
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et81atom.pdb" 
        CNO_to_N = False; S_to_N = False
        folder = ""
        allowed_atoms = ["C","N","O","S"]
    else:
        raise Exception(f"{target} invalid target")


    #### Individual experiment arguments 
    #tag = f"SC{num_unique_supercells}" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
    start_time = -18#-12#-6
    end_time = 18#12#6
    laser_firing_qwargs = dict(
        # pixel sampling method (Neutze) if True - Miller indices if False
        SPI = False,  # sampling method, if False, bragg spots. if True, detector pixels. TODO change name
        SPI_resolution = best_resolution,
        pixels_across = 300,  # for SPI TODO shld go on xfel exp params.
    )
    ##### Crystal params
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 "unit" cells per supercell # Bragg spots will be sampled based on the cell scale, not the supercell scale.
        num_supercells = 1,#100, # 35409
        supercell_simulations = num_unique_supercells, #150
        positional_stdv = 0,#0.2,  #Introduces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        include_symmetries = include_symmetries,  # should unit cell contain symmetries?
        cell_packing = "SC",
        random_waters=NUM_RANDOM_WATER*(RANDOM_WATER==True),
    )
    show_crystal = False

    #### XFEL params
    #TODO make it so reflections don't overwrite same orientation, as stochastic now.
    energy = 7112#7100 # eV
    exp_qwargs = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution), #(best_resolution),#2*np.pi/2
        t_fineness=25,   
        #####crystal stuff (miller)
        max_miller_idx = 25, #None, # = m, [overrides max q so given by q with miller indices (m,m,m)]
        all_miller_indices = True, # False, whether to find all bragg points at or below the max miller index (and between min and max q)
        miller_indices_override=read_hkl(unique_hkl),
        spot_fraction_per_orient=1/cycles_per_bragg_set,
        # first image orientation cardan angles [degrees] 
        num_orients_crys=cycles_per_bragg_set*num_bragg_sets, orientation_axis_crys = [99,99,99],
    )
    same_deviations = False # whether same position deviations between damaged and undamaged crystal 


    # Optional: Choose previous folder for crystal results
    chosen_root_handle = None # None for new. use e.g. "tetra_v1", if want to add images under same params to same results.
    #=========================-------------------------===========================#

    ## DEBUG
    # WARNING we often assume that first crystal is damaged and second is undamaged when plotting. 
    first_crystal_is_damaged = True # True  
    second_crystal_is_damaged = False  # False



    #---------------------------Result handle names---------------------------#
    exp1_qualifier = "real"
    exp2_qualifier = "ideal"
    if chosen_root_handle is None:
        #version_number = 1
        #count = 0
        tag = TAG
        if tag != "":
            tag = "_" + tag        
        #while True:
        # if count > 299:
        #     raise Exception("could not find valid file in " + str(count) + " loops")
        #root_handle = f"{target}{tag}"
        root_handle = f"{target_handle}{tag}"
        exp_name1 = f"{root_handle}_{exp1_qualifier}"
        exp_name2 = f"{root_handle}_{exp2_qualifier}"
        results1_parent_folder = f"{RESULTS_LOCAL_PATH}{exp_name1}/" #_v{version_number}/" 
        results2_parent_folder = f"{RESULTS_LOCAL_PATH}{exp_name2}/" #_v{version_number}/" 
        exp_name1 += f"-{par_idx_for_target}"
        exp_name2 += f"-{par_idx_for_target}"
        
        

        assert(exp_name1!=exp_name2)
        # commented out because parallel
        # if path.exists(path.dirname(results1_parent_folder + exp_name1 + "/")) or path.exists(path.dirname(results2_parent_folder + exp_name2 + "/")):
        #     version_number+=1
        #     count+=1
        #     continue 
        # break
    else:
        exp_name1 = chosen_root_handle + "_" + exp1_qualifier
        exp_name2 = chosen_root_handle + "_" + exp2_qualifier

    #exp_name2 = None

    if SKIP_UNDAMAGED:
        exp_name2 = None # Don't do the undamaged target
    #-------------------------------#
    if SEEDED:
        np.random.seed(0)
        
    src_file_path = inspect.getfile(lambda: None)
    sim_data_dir = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+folder)) + "/"


    # Set up experiments
    #experiment1 = XFEL(exp_name1,energy,**exp_qwargs)
    experiment2 = XFEL(exp_name2,energy,**exp_qwargs)
    # Create Crystals

    crystal = Crystal(pdb_path,allowed_atoms,is_damaged=first_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
    # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
    if same_deviations:
        # we copy the other crystal so that it has the same deviations in coords
        crystal_undmged = copy.deepcopy(crystal)#Crystal(pdb_path,allowed_atoms,cell_dim,is_damaged=False,CNO_to_N = CNO_to_N, **crystal_qwargs)
        crystal_undmged.is_damaged = second_crystal_is_damaged
    else:
        crystal_undmged = Crystal(pdb_path,allowed_atoms,is_damaged=second_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
    if show_crystal:
        crystal.plot_me(300000,water_index = water_index,template="plotly_dark")
    #%
    if laser_firing_qwargs["SPI"]:
        #SPI_result1 = experiment1.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal,results_parent_dir=results1_parent_folder, **laser_firing_qwargs)
        SPI_result2 = experiment2.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal_undmged,results_parent_dir=results2_parent_folder,  **laser_firing_qwargs)
        #stylin(exp_name1,exp_name2,experiment1.max_q,results_parent_dir=results_parent_folder,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2,custom_fig_width=fig_width,custom_fig_height=fig_height)
    else:
        # experiment1.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal, results_parent_dir=results1_parent_folder, **laser_firing_qwargs)
        # exp1_orientations =experiment1.used_orientations
        #rfl_to_sca(exp_name1)
        if exp_name2 != None:
            laser_firing_qwargs["random_orientation"] = False
            #experiment2.set_orientation_set(exp1_orientations)  # pass in orientations to next sim, random_orientation must be false!
            experiment2.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal_undmged, results_parent_dir=results2_parent_folder, **laser_firing_qwargs)
            create_reflection_file(exp_name2,results_parent_dir=results2_parent_folder,artificial_I_scale=1e6/crystal.num_atoms_no_symm())
            _, mtz_file = rfl_to_sca(exp_name2)
            phenix_fcalc(pdb_path,mtz_file)
            phenix_R(pdb_path,mtz_file)

        #stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,results_parent_dir=results_parent_folder, custom_fig_width=fig_width,custom_fig_height=fig_height) # Note we are passing the max q, not max q_scr.


if __name__ == "__main__":
    with Pool(NUM_PARALLEL) as p:
        p.map(main,range(NUM_PARALLEL))
    #all_reflections_to_scalepack() 


# Combine everything into a merged and unmerged scalepack file
if __name__ == "__main__":
    import glob

    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    RESULTS_LOCAL_PATH = "results/"


    num_results = len(target_options)*1#(1+(SKIP_UNDAMAGED==False))
    OldestToLatest = sorted(glob.glob(os.path.join(scattering_dir+RESULTS_LOCAL_PATH, '*/')), key=os.path.getmtime)
    for n in range(num_results):
        all_reflections_to_scalepack(
            OldestToLatest[-n-1].split("/")[-2],
            scattering_dir+RESULTS_LOCAL_PATH,
            out_dir="/home/speno/PhenixWorkspace/data/",
            tag_override="",
            create_mtz=False
        )  

# %%
