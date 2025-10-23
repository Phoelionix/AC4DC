#%%#%%
import matplotlib
#matplotlib.use("pgf")#
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'axes.titlesize':8,     # fontsize of the axes title
    'axes.labelsize':8,    # fontsize of the x and y labels   
    'ytick.labelsize':8,
    'xtick.labelsize':8,
    'legend.fontsize':8,    
    'legend.title_fontsize':8,    
    'lines.linewidth':1,
})
from scatter import *
from multiprocessing import Pool
from sample_scalepack import all_reflections_to_scalepack
import inspect
import datetime
from core_functions import get_sim_params

#TODO auto generate ideal, undamaged, damaged. not undamaged and damaged. (undamaged is called ideal)

NUM_PARALLEL=1
SEEDED = False
RANDOM_WATER = False; NUM_RANDOM_WATER = 0 # 702
DEBUG_WATER = False
QUICK_TEST = False
SKIP_UNDAMAGED = False


target_handle = "copper_sulfate_above_e12_1fs_1" #"copper_sulfate_below_e13_3#"copper_sulfate_above_e12_long_1"#"copper_sulfate_above_e12_14"


def generate_reflections(num_time_points,start_time,end_time): 
    print("=================================================")
    cycles_per_bragg_set = 1 # leave at 1, should be fine
    num_bragg_sets = 1 # 25 # increase this for better stochastic
    num_unique_supercells = 1 # 20
    supercell_scale=1
    num_supercells = 1

    slice_width=0.5 # 99SLICE WIDTH ARBITRARY NOW


    ### Simulate
    #target_options = ["lys_salt","lys_no_salt","neutze","hen","tetra","glycine","fcc"]
    #target_options = ["lys_salt","lys_no_salt"]
    
    
    best_resolution = 0.5 # 1.58 (abdullah) # 2   # resolution (determining max q)
    worst_resolution = 30#30 # 'resolution' corresponding to min q

    #---------------------------------#
    water_index = None # None TODO automate

    include_symmetries=True
    if QUICK_TEST:
        num_bragg_sets = 1
        num_unique_supercells = 1
        include_symmetries=False
    #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_1.4.hkl"

    pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/CuSO4_SC.pdb" 
    CNO_to_N = False; S_to_N = False
    allowed_atoms = ["O","H"]
    folder = ""


    #### Individual experiment arguments 
    #tag = f"SC{num_unique_supercells}" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
    start_time = start_time#-12#-6
    end_time = end_time#12#6
    laser_firing_qwargs = dict(
        # pixel sampling method (Neutze) if True - Miller indices if False
        SPI = False,  # sampling method, if False, bragg spots. if True, detector pixels. TODO change name
        SPI_resolution = best_resolution,
        pixels_across = 300,  # for SPI TODO shld go on xfel exp params.
        do_not_integrate_times=True,
        random_orientation=True,

    )
    ##### Crystal params
    crystal_qwargs = dict(
        supercell_scale = supercell_scale,  # for SC: supercell_scale^3 "unit" cells per supercell # Bragg spots will be sampled based on the cell scale, not the supercell scale.
        num_supercells = num_supercells,#100, # 35409
        supercell_simulations = num_unique_supercells, #150
        positional_stdv = 0,
        zero_bfactors=True,
        #BFACTORS positional_stdv = 0,#0.2,  #Introduces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        include_symmetries = include_symmetries,  # should unit cell contain symmetries?
        cell_packing = "triclinic",
        random_waters=NUM_RANDOM_WATER*(RANDOM_WATER==True),
    )
    show_crystal = False

    #### XFEL params
    #TODO make it so reflections don't overwrite same orientation, as stochastic now.
    energy = 9200#7100 # eV
    exp_qwargs = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution), #(best_resolution),#2*np.pi/2
        t_fineness=num_time_points-1,   
        #####crystal stuff (miller)
        max_miller_idx = 25, #None, # = m, [overrides max q so given by q with miller indices (m,m,m)]
        all_miller_indices = True, # False, whether to find all bragg points at or below the max miller index (and between min and max q)
        spot_fraction_per_orient=None,
        # first image orientation cardan angles [degrees] 
        num_orients_crys=cycles_per_bragg_set*num_bragg_sets, 
        #orientation_axis_crys = [99,99,99],
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
        tag = ""
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
        exp_name1 += f"-slices"
        exp_name2 += f"-slices"
        
        

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
    experiment1 = XFEL(exp_name1,energy,**exp_qwargs)
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
        assert False
    else:
        experiment1.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal, results_parent_dir=results1_parent_folder, **laser_firing_qwargs)
        exp1_orientations =experiment1.used_orientations
        #create_reflection_file(exp_name1,results_parent_dir=results1_parent_folder)
        #rfl_to_sca(exp_name1)
        if exp_name2 != None:
            laser_firing_qwargs["random_orientation"] = False
            experiment2.set_orientation_set(exp1_orientations)  # pass in orientations to next sim, random_orientation must be false!
            experiment2.spooky_laser(start_time,end_time,target_handle,sim_data_dir,crystal_undmged, results_parent_dir=results2_parent_folder, **laser_firing_qwargs)
            #create_reflection_file(exp_name2,results_parent_dir=results2_parent_folder)
            #rfl_to_sca(exp_name2)

    now = datetime.datetime.now().timestamp()
    os.utime(results1_parent_folder[:-1], (now, now))
    if not SKIP_UNDAMAGED:
        os.utime(results2_parent_folder[:-1], (now, now)) # since we are very stupidly and lazily calling scalepack based on most recent modified...
    
        #stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,results_parent_dir=results_parent_folder, custom_fig_width=fig_width,custom_fig_height=fig_height) # Note we are passing the max q, not max q_scr.


if __name__ == "__main__":
    generate_reflections(num_time_points=10,start_time=-1,end_time=1)
    #all_reflections_to_scalepack() 



#%%
import matplotlib
#matplotlib.use("pgf")#
matplotlib.rcParams.update({
    "pgf.texsystem": "pdflatex",
    'font.family': 'serif',
    'text.usetex': True,
    'pgf.rcfonts': False,
    'axes.titlesize':8,     # fontsize of the axes title
    'axes.labelsize':8,    # fontsize of the x and y labels   
    'ytick.labelsize':8,
    'xtick.labelsize':8,
    'legend.fontsize':8,    
    'legend.title_fontsize':8,    
    'lines.linewidth':1,
})
from scatter import *
from multiprocessing import Pool
from sample_scalepack import all_reflections_to_scalepack
import inspect
import datetime
from core_functions import get_sim_params

#TODO auto generate ideal, undamaged, damaged. not undamaged and damaged. (undamaged is called ideal)

NUM_PARALLEL=1
SEEDED = False
RANDOM_WATER = False; NUM_RANDOM_WATER = 0 # 702
DEBUG_WATER = False
QUICK_TEST = False
SKIP_UNDAMAGED = False



# Ideal intensities purely as an extremely inefficient lazy way to scale by pulse profile
def plot_spots_at_time_slices(Intensities, ideal_intensities, miller_indices,times,sqrt=False,normalize=True):
    plt.close()
    stringified = []
    Intensities = copy.deepcopy(Intensities)
    scale_I = np.sum(Intensities[0])
    for elem in miller_indices:
        stringified.append(str(elem[0])+str(elem[1])+str(elem[2])) 
    for I,I_ideal,col, w,  in zip(Intensities,ideal_intensities,['b','g','r'],[1,0.7,0.4]):     
        if normalize:
            I *= scale_I/np.sum(I_ideal) #normalise
        plt.bar(stringified,np.sqrt(I) if sqrt else I,alpha=1,color=col,width=w)
    plt.ylim(0,np.sqrt(np.max(Intensities)))
    plt.xticks(rotation="vertical")
    #plt.show()

def plot_biggest_spots_over_time(Intensities, ideal_intensities, miller_indices,times,ylim=[None,None],sqrt=False,normalize=False,num_points=4,save_path = None,title="",true_normed=False):
    plt.close()
    Intensities = copy.deepcopy(Intensities)
    
  
    scale_I = np.sum(ideal_intensities[0])
    # initialise dict of reflection intensities 
    data = {}
    hkl_keys=[] #ordered keys
    for elem in miller_indices:
        k=str(elem[0])+str(elem[1])+str(elem[2])
        data[k]=[]
        hkl_keys.append(k)
    for i in range(len(Intensities)): # for each time step
        I = Intensities[i]; I_ideal = ideal_intensities[i]
        if normalize:
            I *= scale_I/np.sum(I_ideal)
        if true_normed:
            I/=np.sum(I)
        for hkl, v in zip(hkl_keys,I):
            data[hkl].append(np.sqrt(v) if sqrt else v)
    
        #plt.bar(stringified,np.sqrt(I),alpha=1,color=col,width=w)
    fig_width = 3.49751*0.75 # 20
    fig_height = fig_width*3/4 # 20
    fig,ax = plt.subplots(figsize=(fig_width,fig_height))
    

    Y_list =np.array(list(data.values()))
    sorts = np.argsort(Y_list[:,0])[::-1]
    hkl_list = np.array(list(data.keys()))[sorts][:min(num_points,len(data))]
    Y_list = Y_list[sorts][:min(num_points,len(data))]

    for hkl, Y in zip(hkl_list,Y_list):
        ax.scatter(times,Y,label=hkl,s=4)
        #A = np.vstack([times, np.ones(len(times))]).T
        #m,c =  np.linalg.lstsq()
        
    #ax.set_yscale("log")
    leg =ax.legend(loc="upper left",bbox_to_anchor=[1,1],ncol=1,title="$hkl\\:$")
    leg._legend_box.align = "right"
    ax.set_ylabel("$$I_{hkl}$$")
    if normalize:
        #ax.set_ylabel("$$\\frac{I}{\\Omega}$$")
        ax.set_ylabel("$$I_{hkl}/\\Omega$$")
    if true_normed:
        #ax.set_ylabel("$$\\frac{I}{\\Omega}$$")
        ax.set_ylabel("$$I_{hkl}/\\Sigma I$$")
    assert(not sqrt)
    ax.set_xlabel("Time (fs)")
    ax.set_title(title)
    fig.subplots_adjust(left=0.15, bottom=0.1, right=0.8, top=0.93, wspace=0.3, hspace=0.4) # 3 orb density plots row

    ax.set_ylim(*ylim)
    #plt.ylim(0,None)

    plt.tight_layout(pad=0)

    if save_path is not None:
        plt.savefig(save_path,dpi=300)
    #plt.show()




if __name__ == "__main__":
    import glob
    #target_handle = "copper_sulfate_above_e12_15"
    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    RESULTS_LOCAL_PATH = "results/"


    #num_results = len(times)*(1+(SKIP_UNDAMAGED==False))
    #OldestToLatest = sorted(glob.glob(os.path.join(scattering_dir+RESULTS_LOCAL_PATH, '*/')), key=os.path.getmtime)
    dir = scattering_dir+RESULTS_LOCAL_PATH+target_handle+"_real"
    #print (glob.glob(os.path.join(dir,"*/")))
    #    #for file in glob.glob(dir+"/*"):
    Intensities = []
    ideal_intensities=[]
    files = []
    subdirs =[]
    for filepath in os.listdir(dir):
        subdir = dir+"/"+filepath
        for file in os.listdir(subdir):
            for i, c in  enumerate(filepath):
                if c== "-":
                    break
            files.append(file)
            subdirs.append(subdir)
    assert len(files)==1
    
    for file,subdir in zip(files,subdirs):
        compare_dir=subdir.replace("_real","_ideal")
        dmged_result,undmged_result = get_result(file,subdir,compare_dir=compare_dir)
        Intensities = dmged_result.I
        ideal_intensities = undmged_result.I
        times = dmged_result.T
        miller_indices = dmged_result.miller_indices
    #plot_spots_at_time_slices(Intensities,ideal_intensities,miller_indices,times)
    dir_figures = "../../../output/_Graphs/plots/"
    dir_figures = path.abspath(path.join(__file__ ,dir_figures)) + "/"

    param_dict, param_name_list,unit_list = get_sim_params(target_handle)
    energy,fwhm,photon_count = param_dict["energy"],param_dict["width"],param_dict["fluence"]

    num_points = 10
    note = f"Decay of strongest {num_points} Bragg peaks"
    #title = f"Pulse: "+"$10_power_photon_coun$"+f" ${energy/1000}$"+" keV ph $\\cdot$µm$^"+"{-2}$"+f", ${fwhm}$ fs FWHM"
    title = f"Pulse: "+"$10^{12}$"+f" ${energy/1000}$"+" keV ph $\\cdot$µm$^"+"{-2}$"+f", ${fwhm}$ fs FWHM"
    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,ylim=[2e4,4e4],sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"_high_rfl.png",title=title)
    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,ylim=[0,1.2e4],sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"_low_rfl.png",title=title)
    plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,ylim=[0,None],sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"_rfl.png",title=title)
    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,normalize=True,sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"pulse_scaled_rfl.png",title=title)
    #plot_biggest_spots_over_time(Intensities,ideal_intensities,miller_indices,times,normalize=True,sqrt=False,num_points=num_points,save_path=dir_figures+target_handle+"normed_rfl.png",title=title,true_normed=True)
    assert len(files)>0
    
    # for n in range(num_results):
    #     all_reflections_to_scalepack(
    #         OldestToLatest[-n-1].split("/")[-2],
    #         scattering_dir+RESULTS_LOCAL_PATH,
    #         out_dir="/home/speno/PhenixWorkspace/data/",
    #         tag_override=""
    #     )  


# %%
from core_functions import get_sim_params



# %%
