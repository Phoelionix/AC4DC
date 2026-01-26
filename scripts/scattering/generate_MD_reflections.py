#%%
from scatter import *
from scatter_MD import MD_XFEL,MD_Crystal
from multiprocessing import Pool
from sample_scalepack import all_reflections_to_scalepack
import inspect
import datetime
import sys
sys.path.append('/home/speno/AC4DC')
from scripts.core_functions import get_sim_params
import glob

#TODO Make importable: convert from command line script.

NUM_PARALLEL=1
SEEDED = False
RANDOM_WATER = False; NUM_RANDOM_WATER = 0 # 702
DEBUG_WATER = False
QUICK_TEST = False
SKIP_UNDAMAGED_CONTROL = False; FIRST_UNDAMAGED=False # "second" is the undamaged control # XXX
DEBUG_IGNORE_FE=False
IGNORE_MISSING_SPECIES=True
ELECTRONIC_ONLY=False
NUCLEAR_ONLY=False


assert NUCLEAR_ONLY+ELECTRONIC_ONLY <= 1

PLASMA_SIM_HANDLE_DICT = dict(
    command_line_input=sys.argv[1]
)
MD_results_parent_dir=sys.argv[2]

ground_truth_pdb=None if len(sys.argv)<4 else sys.argv[3]

assert ground_truth_pdb is not None # XXX

target_options = ["command_line_input"]
TAG = "All"
if ELECTRONIC_ONLY:
    TAG="Elect"
if NUCLEAR_ONLY:
    TAG="Nucle"

if FIRST_UNDAMAGED:
    assert SKIP_UNDAMAGED_CONTROL

def generate_single_snapshot_pdb(md_pdb_file,snapshot_index):
    assert False, "Please provide ground truth input - automatic generation unimplemented. "

def main(par_idx,ground_truth_pdb):

    num_time_points = None # use all
    if QUICK_TEST:
        num_time_points=2
    #t_cutoff_frac=0.5
    t_cutoff_frac=None

    fig_width = 3.49751 # 20
    fig_height = fig_width*3/4 # 20
    
    sim_key = target_options[0]
    if par_idx >= NUM_PARALLEL/2 and len(target_options)==2:
        sim_key = target_options[1]
    if NUM_PARALLEL==1:
        par_idx_for_target=1
    else:
        par_idx_for_target = par_idx% int((NUM_PARALLEL/2))+1
    #============------------User params---------==========#
    assert sim_key in target_options
    #sim_key = "lys_no_salt"#"glycine"  #target_options[2]
    #best_resolution = 1.3 # 1.58 (abdullah) # 2   # resolution (determining max q)
    worst_resolution = None#30 # 'resolution' corresponding to min q

    #---------------------------------#
    water_index = None # None TODO automate

    cycles_per_bragg_set = 1
    num_bragg_sets = 1
    num_unique_supercells = 1 
    zero_bfactors=True
    num_supercells=1e5

    if sim_key in PLASMA_SIM_HANDLE_DICT:

        assert num_unique_supercells == 1

        compare_with_undamaged_symmetry=True

        #XXX Put these (unique_hkl + best_resolution, ground_truth_symmetry) as arguments!!!!
        ####XXX#####
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_1.4.hkl"; best_resolution=1.4; ground_truth_symmetry_override="P 43 21 2"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_2.0.hkl"; best_resolution=2.0; ground_truth_symmetry_override="P 43 21 2"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_1.5.hkl"; best_resolution=1.5;ground_truth_symmetry_override="P 21 21 21"
        unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_1.8.hkl"; best_resolution=1.8;ground_truth_symmetry_override="P 21 21 21"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_2.0.hkl"; best_resolution=2.0; ground_truth_symmetry_override="P 21 21 21"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_2.3.hkl"; best_resolution=2.3; ground_truth_symmetry_override="P 21 21 21"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_3.0.hkl"; best_resolution=3.0; ground_truth_symmetry_override="P 21 21 21"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_hemoglobin_6.0.hkl"; best_resolution=6.0; ground_truth_symmetry_override="P 21 21 21"
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_I3C_1.5.hkl"; best_resolution=1.5;ground_truth_symmetry_override=????
        if not compare_with_undamaged_symmetry:
            #ground_truth_pdb="/home/speno/AC4DC/scripts/scattering/targets/4et8H_zero_B.pdb" 
            ground_truth_symmetry_override=None
        ####XXX#####
        
        
        plasma_sim_handle = PLASMA_SIM_HANDLE_DICT[sim_key]
        #pdb_md_snapshots_path = "/home/speno/AC4DC/scripts/scattering/targets/I3C_moldstruct.pdb" 
        #pdb_md_snapshots_path = "/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_moldstruct.pdb" 
        
        # pdb_md_snapshots_path_list = [
        #     "/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_moldstruct1.pdb",
        #     "/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_moldstruct2.pdb",
        #     "/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_moldstruct3.pdb"
        # ]
        # XXX currently just checking for existence of snapshots.pdb to see if results is intended to be used as a sample. Fix this. 
        MD_result_folder_list = [f"{MD_results_parent_dir}/{d}/" for d in os.listdir(MD_results_parent_dir) if os.path.exists(f"{MD_results_parent_dir}/{d}/snapshots.pdb")] # Each file in the directory contains a full MD simulation
        pdb_md_snapshots_path_list=[f"{d}/snapshots.pdb" for d in MD_result_folder_list]
        charges_path_list=[f"{d}/charges.bin" for d in MD_result_folder_list]
        debye_path_list=[f"{d}/debye_data.bin" for d in MD_result_folder_list]
        if ground_truth_pdb is None:
            ground_truth_pdb=generate_single_snapshot_pdb(pdb_md_snapshots_path_list[0],0)
        #ground_truth_pdb="/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_base_structure_moldstruct.pdb" 
        #ground_truth_pdb="/home/speno/AC4DC/scripts/scattering/targets/Lys_salt_moldstruct.pdb" 

        # TODO assert that ground truth has B factors of zero
        CNO_to_N = False; S_to_N = False
        folder = ""
        allowed_atoms = get_sim_elements(plasma_sim_handle)
        if DEBUG_IGNORE_FE and "Fe_fast" in allowed_atoms:
            #allowed_atoms=["N"]
            allowed_atoms.remove("Fe_fast")
        #allowed_atoms = ["C","N","O","S","H"]
        if QUICK_TEST:
            allowed_atoms = ["C"]
    else:
        raise Exception(f"{sim_key} invalid sim_key")


    #### Individual experiment arguments 
    #tag = f"SC{num_unique_supercells}" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + sim_key + tag named according to orientation.
    laser_firing_qwargs = dict(
        # pixel sampling method (Neutze) if True - Miller indices if False
        SPI = False,  # sampling method, if False, bragg spots. if True, detector pixels. TODO change name
        SPI_resolution = best_resolution,
        pixels_across = 300,  # for SPI TODO shld go on xfel exp params.
    )
    ##### Crystal params
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 "unit" cells per supercell # Bragg spots will be sampled based on the cell scale, not the supercell scale.
        num_supercells = num_supercells,#100, # 35409
        supercell_simulations = num_unique_supercells, #150
        #BFACTORS positional_stdv = 0,#0.2,  #Introduces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        cell_packing = "SC",
        random_waters= True if NUM_RANDOM_WATER>0 and RANDOM_WATER else None,
        zero_bfactors=zero_bfactors,
        allow_skip_species=IGNORE_MISSING_SPECIES,
        ignore_water_H=True,
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
        t_fineness=0,   # 1 time point
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
    first_crystal_is_damaged = not FIRST_UNDAMAGED # True  
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
        #root_handle = f"{sim_key}{tag}"
        root_handle = f"{plasma_sim_handle}{tag}"
        exp_name1 = f"{root_handle}_{exp1_qualifier}"
        exp_name2 = f"{root_handle}_{exp2_qualifier}"

        src_file_path = inspect.getfile(lambda: None)
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

    if SKIP_UNDAMAGED_CONTROL:
        exp_name2 = None
    #-------------------------------#
    if SEEDED:
        np.random.seed(0)
        
    src_file_path = inspect.getfile(lambda: None)
    sim_data_dir = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+folder)) + "/"


    # Set up experiments
    experiment1 = MD_XFEL(exp_name1,energy,**exp_qwargs)
    #experiment2 = MD_XFEL(exp_name2,energy,**exp_qwargs)
    experiment2 = XFEL(exp_name2,energy,**exp_qwargs)
    sim_params = get_sim_params(plasma_sim_handle)[0]
    # Create Crystals

    dmged_crystal_targets:list[MD_Crystal] = []
    for snapshots,charges,debye in zip(pdb_md_snapshots_path_list,charges_path_list,debye_path_list): 
        dmged_crystal_targets.append(
            MD_Crystal(num_time_points,snapshots,allowed_atoms,
                       charges,debye,sim_params["start_t"],sim_params["end_t"],
                       t_cutoff_frac=t_cutoff_frac,
                       is_damaged=first_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, 
                       electronic_only=ELECTRONIC_ONLY,nuclear_only=NUCLEAR_ONLY,
                       include_symmetries=False,
                       **crystal_qwargs))
    first_crystal = dmged_crystal_targets[0]

    
    # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
    if same_deviations:
        assert False
        # we copy the other crystal so that it has the same deviations in coords
        crystal_undmged = copy.deepcopy(first_crystal)#Crystal(pdb_md_snapshots_path,allowed_atoms,cell_dim,is_damaged=False,CNO_to_N = CNO_to_N, **crystal_qwargs)
        crystal_undmged.is_damaged = second_crystal_is_damaged
    else:
        # crystal_undmged = MD_Crystal(num_time_points,pdb_md_snapshots_path_list[0],allowed_atoms,
        #                              is_damaged=second_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N,
        #                              **crystal_qwargs)
        crystal_undmged = Crystal(ground_truth_pdb,allowed_atoms,
                                     is_damaged=second_crystal_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N,
                                     include_symmetries=True, 
                                     **crystal_qwargs)
    if show_crystal:
        assert False, "Unimplemented"
        first_crystal.plot_me(300000,water_index = water_index,template="plotly_dark")
    #%
    if laser_firing_qwargs["SPI"]:
        assert False
        #SPI_result1 = experiment1.fire_laser(plasma_sim_handle,sim_data_dir,crystal,results_parent_dir=results1_parent_folder, **laser_firing_qwargs)
        #SPI_result2 = experiment2.fire_laser(plasma_sim_handle,sim_data_dir,crystal_undmged,results_parent_dir=results2_parent_folder,  **laser_firing_qwargs)
        #stylin(exp_name1,exp_name2,experiment1.max_q,results_parent_dir=results_parent_folder,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2,custom_fig_width=fig_width,custom_fig_height=fig_height)
    else:
        I_scale=1e5/crystal_qwargs["num_supercells"]
        experiment1.laser_multi_target(plasma_sim_handle,sim_data_dir,dmged_crystal_targets, results_parent_dir=results1_parent_folder, **laser_firing_qwargs)

        exp1_orientations = experiment1.get_used_orientations()
        create_reflection_file(exp_name1,results_parent_dir=results1_parent_folder,
                               artificial_I_scale=I_scale,symmetry_override=ground_truth_symmetry_override)
        _, mtz_file1 = rfl_to_sca(exp_name1)
        if exp_name2 != None:
            laser_firing_qwargs["random_orientation"] = False
            experiment2.set_orientation_set(exp1_orientations)  # pass in orientations to next sim, random_orientation must be false!
            #experiment2.fire_laser(plasma_sim_handle,sim_data_dir,crystal_undmged, results_parent_dir=results2_parent_folder, **laser_firing_qwargs)
            experiment2.fire_laser(sim_params["start_t"],sim_params["end_t"],plasma_sim_handle,sim_data_dir,crystal_undmged, results_parent_dir=results2_parent_folder, **laser_firing_qwargs)
            create_reflection_file(exp_name2,results_parent_dir=results2_parent_folder,
                                   artificial_I_scale=I_scale,symmetry_override=ground_truth_symmetry_override)
            _, mtz_file2 = rfl_to_sca(exp_name2)
            #fcalc = phenix_fcalc(pdb_md_snapshots_path,best_resolution,real=True)
            phenix_R(ground_truth_pdb,mtz_file2)
        phenix_R(ground_truth_pdb,mtz_file1)
        gen_true_phases=False
        if gen_true_phases: # for e. dens. map making.
            cplx_data = phenix_fcalc(ground_truth_pdb,best_resolution,real=False) 
        fcalc = phenix_fcalc_from_file(ground_truth_pdb,mtz_file1,real=True)
        phenix_R(ground_truth_pdb,fcalc)

    #FIXME
    now = datetime.datetime.now().timestamp()
    os.utime( path.abspath(path.join(__file__ ,"../")) + "/"+ results1_parent_folder[:-1], (now, now))
    if not SKIP_UNDAMAGED_CONTROL:
        os.utime(path.abspath(path.join(__file__ ,"../")) + "/"+ results2_parent_folder[:-1], (now, now)) # since we are very stupidly and lazily calling scalepack based on most recent modified...
    
        #stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,results_parent_dir=results_parent_folder, custom_fig_width=fig_width,custom_fig_height=fig_height) # Note we are passing the max q, not max q_scr.

    print("TODO isomorphous difference map btwn fcalc and ideal")

if __name__ == "__main__":
    if NUM_PARALLEL>1:
        def pooled_func(par_idx):
            main(par_idx,ground_truth_pdb)
        with Pool(NUM_PARALLEL) as p:
            p.map(pooled_func,range(NUM_PARALLEL))
    else:
        main(0,ground_truth_pdb)
    #all_reflections_to_scalepack() 


# Combine everything into a merged and unmerged scalepack file
if __name__ == "__main__":
    import glob

    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    RESULTS_LOCAL_PATH = "results/"
    # TODO FIX 
    ''' 
    num_results = len(target_options)*(1+(SKIP_UNDAMAGED_CONTROL==False))
    OldestToLatest = sorted(glob.glob(os.path.join(scattering_dir+RESULTS_LOCAL_PATH, '*/')), key=os.path.getmtime)
    for n in range(num_results):
        all_reflections_to_scalepack(
            OldestToLatest[-n-1].split("/")[-2],
            scattering_dir+RESULTS_LOCAL_PATH,
            out_dir="/home/speno/PhenixWorkspace/data/",
            tag_override=""
        )  
    '''

# %%
