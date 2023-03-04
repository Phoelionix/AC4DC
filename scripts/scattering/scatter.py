#%%
#TODO make class
#%%
%colors nocolor
import os
os.getcwd() 
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
print(sys.path)

#%%
from Bio.PDB.vectors import Vector as bio_vect
from Bio.PDB.vectors import rotaxis2m
from Bio.PDB.PDBParser import PDBParser
#from sympy.utilities.iterables import multiset_permutations
import itertools
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import cos
from numpy import sin
from plotter_core import Plotter
from scipy.spatial.transform import Rotation as Rotation
import pickle

DEBUG = False 
c_au = 137.036; eV_per_Ha = 27.211385 

class Results():
    def __init__(self,num_points):
        self.azm = np.zeros(num_points)
        self.z = np.zeros(num_points)
        self.q = np.zeros(num_points)
        self.r = np.zeros(num_points)
        self.image_index = np.zeros(num_points)    

    def package_up(self,miller_indices):
        self.q, self.alph = np.meshgrid(self.q,self.azm)  
        self.r, self.alph = np.meshgrid(self.r, self.azm) 
        self.txt = miller_indices 
    def diff(self,other):
        for i in range(len(self.q)):
            subtracted = False
            for j in range(len(other.q)):
                if self.azm[i] == other.azm[j] and self.q[0][i] == other.q[0][j]:
                    self.z[i] = np.abs(self.z[i]- other.z[j])/self.z[i]
                    subtracted = True
            if subtracted == False:
                self.z[i] = -1 
    pass

class Crystal():
    def __init__(self, pdb_fpath, allowed_atoms, rocking_angle = 0.3*np.pi/180, cell_packing = "SC", CNO_to_N = False, orbitals_as_shells = True):
        self.cell_packing = cell_packing
        self.rocking_angle = rocking_angle
        self.sym_factors=[np.array([1,1,1],dtype="float")]
        self.sym_translations = [np.array([0,0,0],dtype="float")]
        self.cell_dim = np.array([1,1,1],dtype="float")              

        parser=PDBParser(PERMISSIVE=1)
        structure_id = os.path.basename(pdb_fpath)
        structure = parser.get_structure(structure_id, pdb_fpath)

        def AC4DC_name(name):  
            """Converts pdb name to names in AC4DC output"""
            if CNO_to_N:
                # light atom approximation. 
                if name == "C" or name == "O":
                    name = "N"
            if  orbitals_as_shells and name not in ["H","He"]:
                name+="_fast"
            if name not in allowed_atoms:
                return None                
            return name
        
        species_dict = {}
        for atom in structure.get_atoms():
            # Get ze data
            R = atom.get_vector()
            name = AC4DC_name(atom.element)
            # Pop into our desired format
            if name == None:
                continue
            if name not in species_dict.keys():
                species_dict[name] = Atomic_Species(name,self) 
            species_dict[name].add_atom(R)
        for str in allowed_atoms:
            if str not in species_dict.keys():
                print("Warning: no",str,"atoms found.")
        self.species_dict = species_dict 

    def set_ff_calculator(self,ff_calculator):
        self.ff_calculator = ff_calculator                  
    def set_cell_dim(self,x,y,z):
        ''' in angstrom'''
        self.cell_dim = np.array([x,y,z])*1.88973
    def add_symmetry(self,symmetry_factor,symmetry_translation):
        self.sym_factors.append(symmetry_factor)
        self.sym_translations.append(symmetry_translation)        

class Atomic_Species():
    def __init__(self,name,crystal):
        self.name = name 
        self.crystal = crystal 
        self.ff = 0 # form_factor
        self.coords = []  # coord of each atom in species in asymmetric unit
    def add_atom(self,vector):
        self.coords.append(vector)
    def set_scalar_form_factor(self,q):
        # F_i(q) = f(q)*T_i(q), where f = self.ff is the time-integrated average 
        self.ff = self.crystal.ff_calculator.f_average(q,self.name) 
    def add_cell_symmetries(self,factor,translation):
        '''Until this is applied, the coords will contain only
        the asymmetric unit. This function adds all symmetries in 
        the unit cell to the coords.            
        '''

        # e.g. (-1,1,-1),(0,0.5,0) -> -X,Y+1/2,-Z  

class XFEL():
    def __init__(self, experiment_name, photon_energy, detector_distance, hemisphere_screen = True, orientation_set = [[0,0,0],], x_orientations = 1, y_orientations = 1,q_cutoff=2.4, pixels_per_ring = 400, num_rings = 50,t_fineness=100):
        """ #### Initialise the imaging experiment's controlled parameters
        experiment_name:
            String that the output folder will be named.        
        photon_energy [eV]:
            Should be the same as that given in the original input file!
        detector_distance [bohr];
            The distance between the target and the centre of the detector 
        q_cutoff [1/bohr]:
            The maximum q allowed.
        pixels_per_ring: 
            The number of different values for phi to plot. phi is the angle that the incident light makes with 
            the y axis, (z-axis is firing line. y-axis is any arbitrary axis perpendicular to z-axis.)
        num_rings: 
            Determines the number of q to calculate. Note pixels are currently just points
        t_fineness:
            Number of time steps to calculate.
        alpha,beta,gamma:
            Angle of rotation about the x, y, and z axis.
        orientation_set:
            contains each set of cardan angles for each orientation to images. Crystal only. TODO replace x_orientations and y_orientations with this.
        """
        self.experiment_name = experiment_name
        self.hemisphere_screen = hemisphere_screen
        self.photon_momentum = 2*np.pi/E_to_lamb(photon_energy)
        self.photon_energy = photon_energy
        self.detector_distance = detector_distance
        self.q_cutoff = q_cutoff
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.x_orientations = x_orientations
        self.y_orientations = y_orientations
        

        self.t_fineness = t_fineness

        self.phi_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)
        self.y_rotation = 0 # Current rotation of crystal (y axis currently)
        self.x_rotation = 0 # Current rotation of crystal (y axis currently)
        self.y_rot_matrix = rotaxis2m(self.y_rotation,bio_vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.x_rotation,bio_vect(1, 0, 0))        
        self.orientation_set = orientation_set 
                      
    def get_ff_calculator(self,start_time,end_time,damage_output_handle):
        ff_calculator = Plotter(damage_output_handle,"y")
        plt.close()
        ff_calculator.initialise_coherence_params(start_time,end_time,self.q_cutoff,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return ff_calculator
    
    def spooky_laser(self, start_time, end_time, damage_output_handle, target, clear_output = False, random_orientation = False, SPI=False):
        """ 
        end_time: The end time of the photon capture in femtoseconds. Not a real thing experimentally, but useful for choosing 
        a level of damage. Explicitly, it is used to determine the upper time limit for the integration of the form factor.
        pdb_fpath: The pdb file's path. Changes the variable self.atoms.
        y_orientations: 
            Number of unique y axis rotations to sample crystal. x_axis_rotations not implemented (yet?).
        random_orientation overrides the XFEL class's orientation_set, replacing each with a random orientation. (get same number of orientations though at present TODO.) 
        """
        ff_calculator = self.get_ff_calculator(start_time,end_time,damage_output_handle)     
        target.set_ff_calculator(ff_calculator)
        self.target = target

        ring = np.empty(self.num_rings,dtype="object")
        
        # Create output folder for results
        directory = "results/" + self.experiment_name + "/"
        print("creating folder:",directory)
        exist_ok = True
        if (os.path.exists(directory)):
            if clear_output:
                exist_ok = False
                for filename in os.listdir(directory):
                    fpath = os.path.join(directory, filename)
                    if os.path.isfile(fpath): 
                        os.remove(fpath)
                os.rmdir(directory)          
        os.makedirs(directory, exist_ok=exist_ok)      

        if SPI:
            result = Results()
            result.z = 0            
            q_samples = np.linspace(0,self.q_cutoff,self.num_rings)
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,bio_vect(1, 0, 0))      
                self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,bio_vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    for i, q in enumerate(q_samples):
                        ring[i] = self.generate_ring(q)
                        #ring[i].I = (ring[i].I+1)
                        #print("q:",q, "x:",ring[i].R,"I[alph=0]",ring[i].I[0])


                    # Initialise stuff that is constant between images (done here to access ring radii.)
                    if rot_y  == 0 and rot_x == 0:
                        azm = self.phi_array
                        radii = np.zeros(self.num_rings)
                        for i in range(len(ring)):
                            radii[i] = ring[i].r           
                        r, alph = np.meshgrid(radii, azm)     
                        q_for_plot, alph = np.meshgrid(q_samples,azm)

                    z = np.zeros(r.shape)  # z is the intensity of the plot colour.
                    for ang in range(len(z)):
                        for pos in range(len(z[ang])):
                            z[ang][pos] = ring[pos].I[ang]                                         
                
                    result.z += z/(self.y_orientations*self.x_orientations)
                self.x_rotation += 2*np.pi/self.x_orientations

        else:
            # Iterate through each orientation of crystal
            used_orientations = []
            for j, cardan_angles in enumerate(self.orientation_set):               
                print("Imaging orientation",j)
                bragg_points, miller_indices,cardan_angles = self.bragg_points(target,cell_packing = target.cell_packing, cardan_angles = cardan_angles,random_orientation=random_orientation)
                used_orientations.append(cardan_angles)
                num_points = int(len(bragg_points))
                result = Results(num_points)
                # Get the q vectors where non-zero
                i = 0
                largest_x = -99999
                for i, G in enumerate(bragg_points):   # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                    if G[0] > largest_x:
                        largest_x = G[0]
                        print("New high G[0]!",G[0])
                    point = self.generate_point(G)
                    #tmp_radii[i] = point.r
                    result.azm[i] = point.phi
                    result.q[i] = point.placeholder_1
                    result.z[i] = point.I
                    result.image_index[i] = j
                result.package_up(miller_indices)
                #Save the result object into its own file within the output folder for the experiment
                fpath = directory + str(cardan_angles) +".pickle"
                pickle_out = open(fpath,"wb")
                pickle.dump(result,pickle_out)
            return used_orientations

    class Feature:
        def __init__(self,q,placeholder_1,theta):
            self.q = q
            self.placeholder_1 = placeholder_1 # magnitude of q parallel to screen
            self.theta = theta          
    class Ring(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    class Spot(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    def generate_ring(self,q):
        '''Returns the intensity(phi) array and the radius for given q.'''
        #print("q=",q)
        #r = self.q_to_r(q)
        #placeholder_1 =self.q_to_q_scr(q)
        theta = self.q_to_theta(q)
        ring = self.Ring(q,placeholder_1,theta)
        ring.I = self.illuminate(ring)
        return ring 
    def generate_point(self,G): # G = vector
        q = np.sqrt(G[0]**2+G[1]**2+G[2]**2)
        #r = self.q_to_r(q)
        placeholder_1 = self.q_to_q_scr(q)# TODO make this r
        if self.hemisphere_screen:
            placeholder_1 = self.q_to_q_scr_curved(G)
        theta = self.q_to_theta(q)
        point = self.Spot(q,placeholder_1,theta)
        point.phi = np.arctan2(G[1],G[0])
        point.G = G

                
        if DEBUG:
            print("placeholder_1",placeholder_1,"phi",point.phi,"G",G[0],G[1],G[2])


        #print("phi",point.phi,"placeholder_1",point.placeholder_1)
        
        phis = np.array([point.phi])
        point.I = self.illuminate(point,phis)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if placeholder_1 != check:
        #     print("Error, placeholder_1 =",placeholder_1,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
    def illuminate(self,feature,phis = None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet."""
        if phis == None:
            phis = self.phi_array
        F = np.zeros(phis.shape,dtype="complex_")
        for species in self.target.species_dict.values():
            species.set_scalar_form_factor(feature.q)
            for R in species.coords:
                    # Rotate to crystal's current orientation 
                    R = R.left_multiply(self.y_rot_matrix)  
                    R = R.left_multiply(self.x_rot_matrix)   
                    # Get spatial factor T
                    T = np.zeros(phis.shape,dtype="complex_")
                    T= self.spatial_factor(phis,R,feature)
                    F += species.ff*T
                    # Rotate atom for next sample            
        I = np.square(np.abs(F))
        

        return I 

    def spatial_factor(self,phi_array,R,feature):
        """ theta = scattering angle relative to z-y plane """ 
        q_z = feature.G[2]#feature.placeholder_1*np.sin(feature.theta) TODO replace with functional method
        q_z = np.tile(q_z,len(phi_array))
        # phi angle of vector relative to x-y plane, pointing from screen centre to point hit. 
        q_y = feature.G[1]#q_z*np.sin(phi_array)
        q_x = feature.G[0]#q_z*np.cos(phi_array)
        q_vect = np.column_stack([q_x,q_y,q_z])
        for i in range(len(self.target.sym_factors)):
            coord = R.get_array()
            coord = np.multiply(R.get_array(),self.target.sym_factors[i]) + np.multiply(self.target.cell_dim,self.target.sym_translations[i])
            spatial_structure_factor = np.exp(-1j*np.dot(q_vect,coord))   
        return spatial_structure_factor
    def bragg_points(self,crystal, cell_packing, cardan_angles,random_orientation=False):
        ''' 
        Using the unit cell structure, find non-zero values of q for which bragg 
        points appear.
        lattice_vectors e.g. = [a,b,c] - length of each spatial vector in orthogonal basis.
        '''

        def get_G(miller_indices,cartesian=True):
            '''
            Get G in global cartesian coordinates (crystal axes not implemented),
            where G = G[0]x + G[1]y + G[2]z
            '''       #TODO vectorise somehow      
            G = np.zeros(miller_indices.shape)
            if cartesian:
                for i in range(len(miller_indices)):
                    G[i] = np.array(np.dot(miller_indices[i],b))
                    if (i == 2 or i == 10) and DEBUG:
                        print("G Debug")
                        print(b)
                        print(miller_indices[i],G[i])
                #G = G.reshape(3,)
                G, used_angles = self.rotate_G_to_orientation(G,*cardan_angles,random = random_orientation) 
            else: 
                G = "non-cartesian Unimplemented"
            return G, used_angles
        
        # (h,k,l) is G (subject to selection condition) in lattice vector basis (lattice vector length = 1 in each dimension):
        # a = lattice vectors, b = reciprocal lattice vector
        if cell_packing == "SC" or cell_packing == "FCC" or cell_packing == "BCC" or cell_packing == "FCC-D":
            if cell_packing == "SC":
                a = np.array([[1,0,0],[0,1,0],[0,0,1]])
            if cell_packing == "BCC":
                a = 0.5*np.array([[-1,1,1],[1,-1,1],[1,1,-1]])
            if cell_packing == "FCC":
                a = 0.5*np.array([[0,1,1],[1,0,1],[1,1,0]])                
            a = np.multiply(a,crystal.cell_dim) # idk if this is correct for rhomboids TODO
        b1 = np.cross(a[1],a[2])
        b2 = np.cross(a[2],a[0])
        b3 = np.cross(a[0],a[1])
        b = 2*np.pi*np.array([b1,b2,b3])/(np.dot(a[0],np.cross(a[1],a[2])))
        
        # Cast a wide net, catching all possible permutations of miller indices.
        #   G = hb1 + kb2 + lb3. (bi = lattice vector)
        #   !Attention! Assuming vectors are orthogonal.
        q_1_max = self.q_cutoff
        q_2_max = self.q_cutoff
        q_3_max = self.q_cutoff        # q = (0,0,l) case.
        h_max = 0
        while h_max*np.sqrt(sum(pow(element, 2) for element in b[0])) <= q_1_max:
            h_max += 1
        k_max = 0
        while k_max*np.sqrt(sum(pow(element, 2) for element in b[1])) <= q_2_max:
            k_max += 1 
        l_max = 0
        while l_max*np.sqrt(sum(pow(element, 2) for element in b[2])) <= q_3_max:
            l_max += 1                  
        h_set = np.arange(-h_max,h_max+1,1)
        k_set = np.arange(-k_max,k_max+1,1)
        l_set = np.arange(-l_max,l_max+1,1)
        indices = list(itertools.product(h_set,k_set,l_set))
        indices = np.array([*set(indices)])   # * expands the set into the list. Faster than using List() apparently.
    
        # Selection rules/Fish bait
        if cell_packing == "SC":
            selection_rule = lambda f: True
        if cell_packing == "BCC":
            selection_rule = lambda f: (f[0]+f[1]+f[2])%2==0  # All even
        if cell_packing == "FCC":
            selection_rule = lambda f: (np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) in [0,3]   # All odd or all even.
        if cell_packing == "FCC-D":
            selection_rule = lambda f: (np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) == 3 or ((np.abs(f[0])%2+np.abs(f[1])%2+np.abs(f[2])%2) == 0 and (f[0] + f[1] + f[2])%4 == 0)  # All odd or all even.    

        # Set G, and retrieve the cardan angles used (in case of random orientations)
        G_temp, cardan_angles = get_G(indices) 
        random_orientation = False # Only do random orientations once

        # Purely a matter of truncation
        q_cutoff_rule = lambda g: np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2))<= self.q_cutoff
        #q_cutoff_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.q_cutoff
        q0 = self.photon_momentum
        # Mosaicity:
        
        # Catch the miller indices with a boolean mask
        mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(q_cutoff_rule,1,G_temp)*np.apply_along_axis(self.mosaicity,1,G_temp)
        indices = indices[mask]

        print("Cardan angles:",cardan_angles)
        print("Number of points:", len(indices))   
        for elem in indices:
            if DEBUG:
                print(elem)
        

        G, cardan_angles = get_G(indices,cardan_angles)
        if DEBUG:
            print(G)
        return G, indices, cardan_angles

    def rotate_G_to_orientation(self,G,alpha=0,beta=0,gamma=0,random=False):
        '''
        Rotates the G vectors to correspond to the crystal when oriented to the given cardan angles. 
        alpha: x-axis angle [radians].
        beta:  y-axis angle [radians].
        gamma: z-axis angle [radians].
        '''
        # Handle case of shape = (N,3)
        transposed = False
        if G.shape[0] != 3:
            G = G.T
            transposed = True

        if random == True:
            alpha = np.random.random_sample()*2*np.pi
            beta = np.random.random_sample()*2*np.pi
            gamma = np.random.random_sample()*2*np.pi

        R = Rotation.from_euler('xyz',[alpha,beta,gamma])
        # Apply rotation
        print("Rotation matrix:")
        print(R.as_matrix())
        G = np.matmul(np.around(R.as_matrix(),decimals=10),G)

        if transposed:
            G = G.T
        return G, [alpha,beta,gamma]
                
    # q_abs [1/bohr]
    def q_to_theta(self,q):
        lamb = E_to_lamb(self.photon_energy)
        theta = np.arcsin(lamb*q/(4*np.pi)) 
        if theta < 0 or theta > 45:
            print("warning, theta out of bounds")        
        return theta


    # def q_to_r(self, q):
    #     """Returns the distance from centre of screen that is struck by photon which transferred q (magnitude)"""
    #     D = self.detector_distance #bohr
    #     theta = self.q_to_theta(q)
    #     return D*np.tan(2*theta)

    def q_to_scr_pos(self,q):
        """ 
        """
        theta = self.q_to_theta(q)
        # From here it will help to consider q as referring to the final momentum.
        v_x = c_au*np.cos(2*theta)
        v_y = c_au*np.sin(2*theta)
        D = self.detector_distance
        t = D/v_x
        radius = v_y*t
        return radius
    def q_to_q_scr(self,q):
            """ Returns the screen-parallel component of q (or equivalently the final momentum).
            """
            theta = self.q_to_theta(q)
            return q*np.cos(theta)
    def q_to_q_scr_curved(self,G):
        return np.sqrt(G[0]**2+G[1]**2)
    
    def mosaicity(self,q_vect):
        ''' determines whether vector q is allowed.'''
        # Rocking angle version
        q = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + q_vect[2]**2)
        theta = self.q_to_theta(q)
        k = self.photon_momentum
        #Assuming theta > 0:
        delt = self.target.rocking_angle
        min_err = k*sin(theta - delt/2)
        max_err = k*sin(theta+delt/2)
        err = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + (q_vect[2]-k)**2) - k
        if min_err <= err <= max_err:
            return True
        return False
        

    def r_to_q(self,r):
        D = self.detector_distance
        lamb = E_to_lamb(self.photon_energy) 
        theta = np.arctan(r/D)/2
        q = 4*np.pi/lamb*np.sin(theta)
        return q 
    # def r_to_q_scr(self,r):
    #     q = self.r_to_q(r)
    #     return self.q_to_q_scr(q)

def E_to_lamb(photon_energy):
    """Energy to wavelength in A.U."""
    E = photon_energy  # eV
    return 2*np.pi*c_au/(E/eV_per_Ha)
    
    # Crystalline.

        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.

def load_results(directory):
    for filename in os.listdir(directory):
        fpath = os.path.join(directory, filename)
        if os.path.isfile(fpath):
            pickle.load(fpath)

def scatter_plot(result_handle, compare_handle = None, cmap_power = 1, cmap = None, min_alpha = 0.05, max_alpha = 1, solid_colour = "white", show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False,  fixed_dot_size = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None, **cmesh_kwargs):
    '''
    Plots the simulated scattering image.
    result_handle:

    compare_handle:
        results/compare_handle/ is the directory of results that will be subtracted from those in the results directory.
        Must have same orientations.
    '''
    results_dir = "results/"+result_handle+"/"
    compare_dir = None
    if compare_handle!= None:
        compare_dir = "results/"+compare_handle+"/"
    kw = cmesh_kwargs
    # if "color" not in cmesh_kwargs: 
    #     kw["color"] = 'k'
    # if "ls" not in cmesh_kwargs.keys():
    #     kw["ls"] = 'none'
    fig = plt.figure()
      
    #ax = Axes3D(fig) 
    ax = fig.add_subplot(projection="polar")
    # Iterate through each file to get normalisation values
    max_z = -np.inf
    min_z = np.inf
    for filename in os.listdir(results_dir):
        fpath = os.path.join(results_dir, filename)
        if os.path.isfile(fpath):
            result = pickle.load(open(fpath,'rb'))
        if compare_dir != None:
            if filename in os.listdir(compare_dir):
                fpath2 = os.path.join(compare_dir, filename)
                result2 = pickle.load(open(fpath2,'rb'))
                result.diff(result2)
            else:
                break
        max_z = max(max_z,np.max(result.z))
        min_z = min(min_z,np.min(result.z))
    if log_I:
        max_z = np.log(max_z)
        min_z = np.log(min_z)
    # Iterate through each orientation (one for each file) 
    for filename in os.listdir(results_dir):
        fpath = os.path.join(results_dir, filename)
        if os.path.isfile(fpath):
            result = pickle.load(open(fpath,'rb'))
        if compare_dir != None:
            if filename in os.listdir(compare_dir):
                fpath2 = os.path.join(compare_dir, filename)
                result2 = pickle.load(open(fpath2,'rb'))
                result.diff(result2)
            else:
                # Not bothering to check if all orientations of result2 is in result 1 
                print("ERROR, missing matching orientation in second results directory")
                break            

        radial_axis = result.r
        if plot_against_q:
            radial_axis =  result.q
        if len(result.z.shape) == 1:
            radial_axis = radial_axis[0]
            ## Point-like (Crystalline)
            # get points at same position on screen.
            identical_count = np.zeros(result.z.shape)
            tmp_z = np.zeros(result.z.shape)
            processed_copies = []
            unique_values_mask = np.zeros(result.z.shape,dtype="bool")
            for i in range(len(result.z)):
                if (radial_axis[i], result.azm[i],result.image_index[i])  in processed_copies:
                    unique_values_mask[i] = False
                    continue
                else:
                    unique_values_mask[i] = True
                    
                # Boolean masks
                matching_rad_idx = np.nonzero(radial_axis == radial_axis[i])  # numpy note: equiv. to np.where(condition). Non-zero part irrelevant.
                matching_alph_idx = np.nonzero(result.azm == result.azm[i])
                for elem in matching_rad_idx[0]:
                    if elem in matching_alph_idx[0]:
                        identical_count[i] += 1
                        matching = result.z[(radial_axis == radial_axis[i])*(result.azm == result.azm[i])]
                        for value in matching:
                            tmp_z[i] += value
                #processed_copies.append((radial_axis[i],result.azm[i],result.image_index[i]))
            #print("processed copies:",processed_copies)
            
            identical_count = identical_count[unique_values_mask]
            tmp_z = tmp_z[unique_values_mask]
            radial_axis = radial_axis[unique_values_mask]
            azm = result.azm[unique_values_mask]
            
            z = tmp_z

            if log_I: 
                z = np.log(z)
            else:
                z = z           

            #debug_mask = (0.01 < radial_axis[0])*(radial_axis[0] < 100)
            #print(identical_count[debug_mask])
            #print(radial_axis[0][debug_mask])
            #print(azm[debug_mask ]*180/np.pi)

            import matplotlib as mpl
            from matplotlib import cm

            #https://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
            
            from matplotlib.colors import to_rgb
            if cmap != None:
                class MplColorHelper:

                    def __init__(self, cmap_name, start_val, stop_val):
                        self.cmap_name = cmap_name
                        self.cmap = plt.get_cmap(cmap_name)
                        self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
                        self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

                    def get_rgb(self, val):
                        val = ((val-min_z)/(max_z-min_z))**cmap_power
                        if val < 0 or val > 1:
                            print("error in val!")
                        return self.scalarMap.to_rgba(val)
                COL = MplColorHelper(cmap, 0, 1) 
                thing = np.array([[COL.get_rgb(K)[0], COL.get_rgb(K)[1],COL.get_rgb(K)[2],np.clip((K*(max_alpha-min_alpha))/(max_z) + min_alpha,min_alpha,None)] for K in z])
                #print("THING",thing)
                colours = [(r,g,b,a) for r,g,b,a in thing] 
            else:
                r,g,b = to_rgb(solid_colour)
                colours = [(r,g,b,a) for a in np.clip(z/max_z,min_alpha,max_alpha)]
            alpha = None

            # Dot size
            dot_param = z
            if crystal_pattern_only:
                dot_param = identical_count
            if not log_dot:
                dot_param = np.e**(dot_param)     
            norm = np.max(dot_param)
            s = [100*dot_size*x/norm for x in dot_param]
            if fixed_dot_size:
                s = [100*dot_size for x in dot_param]
            ax.scatter(azm,radial_axis,c=colours,s=s,alpha=alpha,**kw)
            plt.grid(False) 
            #TODO only show first index or something. or have option to switch between image indexes. oof coding.
            if show_labels:
                for i, txt in enumerate(result.txt[unique_values_mask]):
                    ax.annotate("%.0f" % txt[0]+","+"%.0f" % txt[1]+","+"%.0f" % txt[2]+",", (azm[i], radial_axis[i]),ha='center')        
        else:
            ## Continuous (SPI)
            if log_I: 
                z = np.log(result.z)
            else:
                z = result.z        
            if log_I and cutoff_log_intensity != None:
                #cutoff_log_intensity = -1
                z -= cutoff_log_intensity
                z[z<0] = 0
                pass
            ax.pcolormesh(result.alph, radial_axis, z,**kw)
            ax.plot(result.azm, radial_axis, color = 'k',ls='none')
            plt.grid()  # Make the grid lines represent one unit cell (when implemented).

        if radial_lim:
            bottom,top = plt.ylim()
            plt.ylim(bottom,radial_lim)
        if log_radial:
            plt.yscale("log")  
        ax.set_facecolor("black")

#%% 
#  #TODO make this nicer and pop in a function 
DEBUG = False
#energy = 6000 # Tetrapeptide 
energy =17445   #crambin - from resolutions. in pdb file, need to double check calcs: #q_min=0.11,q_cutoff = 3.9,  my defaults: pixels_per_ring = 400, num_rings = 200

exp_name1 = "Lys_9.95_random"
exp_name2 = "lys_9.80_random"
#orientation_set = [[np.pi/4,0,0],[0,0,0],[0,np.pi/4,0]]
#orientation_set = [[0,0,z_angle] for z_angle in np.linspace(0,2*np.pi,20,endpoint=False)]
orientation_set = [[0,0,0] for g in range(25)]
#orientation_set = [[6.153223907916695, 3.9691232140013524, 2.2006805479640956]]
experiment1 = XFEL(exp_name1,energy,100, hemisphere_screen = False, orientation_set = orientation_set, q_cutoff=10, pixels_per_ring = 400, num_rings = 400,t_fineness=100)
experiment2 = XFEL(exp_name2,energy,100, hemisphere_screen = False, orientation_set = orientation_set, q_cutoff=10, pixels_per_ring = 400, num_rings = 400,t_fineness=100)

random_orientation = True


sym_translations = [np.array([0,0,0])]
cell_dim = [np.array([1,1,1])]  

pdb_path = "/home/speno/AC4DC/scripts/scattering/3u7t.pdb" #Crambin

allowed_atoms_1 = ["C_fast","N_fast","O_fast","S_fast"]
allowed_atoms_2 = ["C_fast","N_fast","O_fast","S_fast"]
end_time_1 = -9.95
end_time_2 = -9.80
output_handle = "Improved_Lys_mid_6"

crystal = Crystal(pdb_path,allowed_atoms_1,CNO_to_N=False,cell_packing = "SC")
crystal.set_cell_dim(20, 20, 20)
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0]))

SPI = False

exp1_orientations = experiment1.spooky_laser(-10,end_time_1,output_handle,crystal, random_orientation = random_orientation, SPI=SPI)

experiment2.orientation_set = exp1_orientations
experiment2.spooky_laser(-10,end_time_2,output_handle,crystal, random_orientation = False, SPI=SPI)

#%%
# stylin' 
experiment1_name = exp_name1
experiment2_name = exp_name2 
#experiment2_name = None

#####


font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 24}

plt.rc('font', **font)

use_q = True
log_radial = False
log_I = True
cutoff_log_intensity = -1#-1
import colorcet as cc; import cmasher as cmr
cmap = "plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
cmap_power = 1.6
min_alpha = 0.1
max_alpha = 1
colour = "y"

# screen_radius = 150#55#165    #
# q_scr_lim = experiment.r_to_q(screen_radius) #experiment.r_to_q_scr(screen_radius)#3.9  #NOTE this won't be the actual max placeholder_1 but ah well.
# zoom_to_fit = True
# ####### n'
# # plottin'

# if zoom_to_fit:
#     radial_lim = screen_radius#min(screen_radius,experiment.q_to_r(experiment.q_cutoff))
#     print(radial_lim)
#     if use_q:
#         #radial_lim = experiment.r_to_q_scr(radial_lim)
#         radial_lim = q_scr_lim #radial_lim = min(q_scr_lim,experiment.q_to_q_scr(experiment.q_cutoff))
#         print(radial_lim)
# else:
#     radial_lim = None

radial_lim = 5#10 #TODO fix above then remove this

scatter_plot(experiment1_name, experiment2_name, fixed_dot_size = True, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
fig = plt.gcf()
fig.set_figwidth(20)
fig.set_figheight(20)

#NEED TO CHECK. We have a 1:1 mapping from q to q_parr, but with our miller indices we are generating multiple q_parr with diff q.
# So we SHOULD get the same q_parr with different q_z. Which makes sense since we are just doing cosine. But still icky maybe?
# Need to double check we get different intensities for same q_parr. Pretty sure that's implemented.


#TODO 
# log doesnt work atm.
# %%


# get sinc functions with non-infinite crystal
# 