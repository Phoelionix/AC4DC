from scatter import res_to_q

print("====================")
print("Refreshed imaging params")
print("====================")

top_resolution = 2
bottom_resolution = None
default_dict = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = False,
        pixels_across = 100,  # for SPI, shld go on xfel params.
        random_orientation = True,  # NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        cell_scale = 3,  # for SC: cell_scale^3 unit cells 
        positional_stdv = 0.2, # RMS in atomic coord position [angstrom] (set to 0 below if crystal, since rocking angle handles this aspect)
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 1.2,  # (approximating mosaicity)
        orbitals_as_shells = True,
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
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
        #crystallographic orientations (not consistent with SPI yet)
        # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations      
        num_orients_crys=15,
        orientation_axis_crys = None,#[1,1,0]        
    ),
)
lys_dict = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = False,
        pixels_across = 100,  # for SPI, shld go on xfel params.
        random_orientation = True,  # NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        cell_scale = 3,  # for SC: cell_scale^3 unit cells 
        positional_stdv = 0.2, # RMS in atomic coord position [angstrom] (set to 0 below if crystal, since rocking angle handles this aspect)
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 0.15,  # (approximating mosaicity)
        orbitals_as_shells = True,
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
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
        #crystallographic orientations (not consistent with SPI yet)
        # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations      
        num_orients_crys=20,
        orientation_axis_crys = None,#[1,1,0]        
    ),
)
SPI_default_dict = dict( 

    #### Individual experiment arguments
    laser = dict(
        SPI = True,
        pixels_across = 100,  # for SPI, shld go on xfel params.
        random_orientation = True,  # NB: orientation is synced with undamaged crystal imaging  (TODO random orientation should not be separate from the other orient params...)
    ),
    ##### Crystal params
    crystal = dict(
        cell_scale = 1,  # for SC: cell_scale^3 unit cells 
        positional_stdv = 0.2, # RMS in atomic coord position [angstrom] (set to 0 below if crystal, since rocking angle handles this aspect)
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        rocking_angle = 0.1,  # (approximating mosaicity)
        orbitals_as_shells = True,
        #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
    ),
    #### XFEL params    
    experiment = dict(
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
        #crystallographic orientations (not consistent with SPI yet)
        # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Incompatible with random_orientations      
        num_orients_crys=3,
        orientation_axis_crys = None,#[1,1,0]        
    ),
)