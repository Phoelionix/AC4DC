import scatter 

### Simulate
    #TODO
    # implement rhombic miller indices as the angle is actually 120 degrees on one unit cell lattice vector (or just do SPI)
def main():
    target_options = ["neutze","hen","tetra"]
    #============------------User params---------==========#

    target = "tetra" #target_options[2]
    top_resolution = 2
    bottom_resolution = None#30

    #### Individual experiment arguments 
    start_time = -6
    end_time = 6#10 #0#-9.80    
    laser_firing_qwargs = dict(
        SPI = True,
        pixels_across = 100,  # for SPI, shld go on xfel params.
    )
    ##### Crystal params
    crystal_qwargs = dict(
        cell_scale = 2,  # for SC: cell_scale^3 unit cells 
        positional_stdv = 0.2, # RMS in atomic coord position [angstrom] (set to 0 below if crystal, since rocking angle handles this aspect)
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 0.1,  # (approximating mosaicity)
        orbitals_as_shells = True,
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    )

    #### XFEL params
    tag = "" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
    #TODO
    random_orientation = True # if true, overrides orientation_set with random orientations (only affects non-SPI)
    energy = 12000 # eV
    exp_qwargs = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(bottom_resolution),#None #angstrom
        q_cutoff = res_to_q(top_resolution),#2*np.pi/2
        t_fineness=100,   
        #####crystal stuff
        max_triple_miller_idx = None, # = m, where max momentum given by q with miller indices (m,m,m)
        ####SPI stuff
        num_rings = 20,
        pixels_per_ring = 20,
        # first image orientation (cardan angles, degrees) 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        ######
    )

    #crystallographic orientations (not implemented for SPI yet)
    num_orients = 15
    # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Overridden if random orientations.
    ax_x = 1
    ax_y = 1
    ax_z = 0


    # Optional: Choose previous folder for crystal results
    chosen_root_handle = None # None for new. use e.g. "tetra_v1", if want to add images under same params to same results.
    #=========================-------------------------===========================#

    #----------------------- Turn off stdv for crystal -----------------------#
    if laser_firing_qwargs["SPI"] == False:
        crystal_qwargs["positional_stdv"] = 0
    #---------------------------Result handle names---------------------------#
    exp1_qualifier = "_real"
    exp2_qualifier = "_ideal"
    if chosen_root_handle is None:
        version_number = 1
        count = 0
        while True:
            if tag != "":
                tag = "_" + tag
            if count > 99:
                raise Exception("could not find valid file in " + str(count) + " loops")
            results_dir = "results/" # needs to be synced with other functions
            root_handle = str(target) + tag + "_v" + str(version_number)
            exp_name1 = root_handle + "_" + exp1_qualifier + "real"
            exp_name2 = root_handle + "_" + exp2_qualifier + "ideal"    
            if os.path.exists(os.path.dirname(results_dir + exp_name1 + "/")) or os.path.exists(os.path.dirname(results_dir + exp_name2 + "/")):
                version_number+=1
                count+=1
                continue 
            break
    else:
        exp_name1 = chosen_root_handle + "_" + exp1_qualifier + "real"
        exp_name2 = chosen_root_handle + "_" + exp2_qualifier + "ideal"  

    #exp_name2 = None
    #----------orientation thing that needs to be to XFEL class (TODO)-----------------------#

    if ax_x == ax_y == ax_z and ax_z == 0 and random_orientation == False:
        raise Exception("axis vector has 0 length")
    # if not random_orientation:
    orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
    orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
    orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
    print("ori set if not random:", orientation_set)  #TODO double check shows when random and if so disable when random
    #---------------------------------#
    if target == "neutze": #T4 virus lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/2lzm.pdb"
        cell_dim = np.array((61.200 ,  61.200 ,  61.2)) #TODO should read this automatically...
        target_handle = "D_lys_neutze_simple_7"  
        folder = ""
        allowed_atoms_1 = ["N_fast","S_fast"]
        CNO_to_N = True
    elif target == "hen": # egg white lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb"
        cell_dim = np.array((79.000  , 79.000  , 38.000))
        target_handle = "D_lys_neutze_simple_7"  
        folder = ""
        allowed_atoms_1 = ["N_fast","S_fast"]
        CNO_to_N = True
    elif target == "tetra": 
        pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb" 
        cell_dim = np.array((4.813 ,  17.151 ,  29.564))
        target_handle = "12-5-2_tetra_2" #"carbon_12keV_1" # "carbon_6keV_1" #"carbon_gauss_32"
        folder = "tetra"
        allowed_atoms_1 = ["N_fast"]
        CNO_to_N = True
        #CNO_to_N = False
    else:
        raise Exception("'target' invalid")
    #-------------------------------#

    import inspect
    src_file_path = inspect.getfile(lambda: None)
    parent_directory = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+folder)) + "/"


    # Set up experiments
    experiment1 = XFEL(exp_name1,energy,orientation_set = orientation_set,**exp_qwargs)
    experiment2 = XFEL(exp_name2,energy,orientation_set = orientation_set,**exp_qwargs)
    # Create Crystals

    crystal = Crystal(pdb_path,allowed_atoms_1,cell_dim,is_damaged=True,CNO_to_N = CNO_to_N, **crystal_qwargs)
    # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
    # we copy the other crystal so that it has the same deviations in coords
    crystal_undmged = copy.deepcopy(crystal)#Crystal(pdb_path,allowed_atoms_1,cell_dim,is_damaged=False,CNO_to_N = CNO_to_N, **crystal_qwargs)
    crystal_undmged.is_damaged = False
    crystal.plot_me(20000)

    if laser_firing_qwargs["SPI"]:
        SPI_result1 = experiment1.spooky_laser(start_time,end_time,target_handle,parent_directory,crystal, random_orientation = random_orientation, **laser_firing_qwargs)
        SPI_result2 = experiment2.spooky_laser(start_time,end_time,target_handle,parent_directory,crystal_undmged, random_orientation = False, **laser_firing_qwargs)
        stylin(SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
    else:
        exp1_orientations = experiment1.spooky_laser(start_time,end_time,target_handle,parent_directory,crystal, random_orientation = random_orientation, **laser_firing_qwargs)
        create_reflection_file(exp_name1)
        if exp_name2 != None:
            experiment2.orientation_set = exp1_orientations  # pass in orientations to next sim, random_orientation must be false so not overridden!
            experiment2.spooky_laser(start_time,end_time,target_handle,parent_directory,crystal_undmged, random_orientation = False, **laser_firing_qwargs)
        stylin()


