from scatter import res_to_q
import copy

print("====================")
print("Refreshed imaging params")
print("====================")

best_resolution = 2
worst_resolution = None
default_dict = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = False,
        SPI_resolution = best_resolution,
        pixels_across = 50,  # for SPI, shld go on xfel params.
        random_orientation = True,  # infinite cryst sim only, NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        cell_scale = 2,  # for SC: cell_scale^3 unit cells 
        positional_stdv = 0, # Introduces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out.
        include_symmetries = False,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution),#2*np.pi/d
        t_fineness=25,   
        #####crystal stuff
        max_triple_miller_idx = None, # = m, where max momentum given by q with miller indices (m,m,m)
        ####SPI stuff
        num_rings = 20,
        pixels_per_ring = 25,
        # first image orientation (cardan angles, degrees) 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        #crystallographic orientations (not consistent with SPI yet)
        num_orients_crys=20,
        orientation_axis_crys = None,#[1,1,0]   # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations   
    ),
)
default_dict_SPI = copy.deepcopy(default_dict)
default_dict_SPI["laser"]["SPI"] = True 
# Water background
background_dict = copy.deepcopy(default_dict_SPI)
#background_dict["crystal"]["positional_stdv"] = 10  ### Idea: When we average over the distributions, this will form a background. 
background_dict["crystal"]["positional_stdv"] = 0  ### Idea:  generate multiple solvents
background_dict["crystal"]["cell_scale"] = 1
background_dict["crystal"]["include_symmetries"] = False 
background_dict["experiment"]["t_fineness"] = 10  # We don't need to be as accurate with this.
##