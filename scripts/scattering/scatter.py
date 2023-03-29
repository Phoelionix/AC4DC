'''
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================
(C) Spencer Passmore 2023
'''
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
from Bio.PDB.vectors import Vector as Bio_Vect
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
from matplotlib.colors import to_rgb
import matplotlib as mpl
from matplotlib import cm
import copy
import pickle

DEBUG = False 
c_au = 137.036; eV_per_Ha = 27.211385 

class Results():
    def __init__(self,num_points):
        self.phi = np.zeros(num_points)
        self.phi_aligned = np.zeros(num_points)
        self.I = np.zeros(num_points)
        self.q = np.zeros(num_points)
        self.r = np.zeros(num_points)
        self.image_index = np.zeros(num_points)    

    def package_up(self,miller_indices):
        _, self.phi_mesh = np.meshgrid(self.q,self.phi)  
        self.r, self.phi_mesh = np.meshgrid(self.r, self.phi)
        self.q, self.phi_aligned_mesh = np.meshgrid(self.q,self.phi_aligned) 
        self.miller_indices = miller_indices 
    def diff(self,other):
        '''
        Get single-point R factors between two images
        '''
        self.R = np.zeros(self.I.shape)
        for i in range(len(self.q)):
            subtracted = False
            for j in range(len(other.q)):
                if self.phi[i] == other.phi[j] and self.q[0][i] == other.q[0][j]:   # Using phi is equivalent to using phi_aligned, as this function is only used for same-orientation comparisons.
                    #self.I[i] = np.abs(self.I[i]- other.I[j])/self.I[i] # Intensity
                    self.R[i] = np.abs(np.sqrt(self.I[i])- np.sqrt(other.I[j]))/np.sqrt(max(self.I[i],other.I[j]))# form factor  # TODO not sure why later times is sometimes larger but probs interference thing. Hopefully will disappear when we use large time gaps... or at least when we average over to get R factor.
                    subtracted = True
            if subtracted == False:
                self.R[i] = -1 
class Results_SPI():
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
        
        ## Dictionary for going from pdb to ac4dc names.
        # different names
        PDB_to_AC4DC_dict = dict(
            NA = "Sodion", CL = "Chloride",
        )
        # Same names
        for elem in ["H","He","C","N","O","P","S"]:
            PDB_to_AC4DC_dict[elem] = elem
        # Modify the dictionary according to arguments.
        for k,v in PDB_to_AC4DC_dict.items():
            # Light atom approximation. 
            if CNO_to_N:
                if v in ["C","O"]: 
                    v = "N"                    
            # Orbital approximation.
            if orbitals_as_shells:
                if v not in ["H","He","Sodion","Chloride"]:
                    v = v+"_fast"      
            PDB_to_AC4DC_dict[k] = v  
        species_dict = {}
        pdb_atoms = []
        pdb_atoms_ignored = ""
        for atom in structure.get_atoms():
            # Get ze data
            R = atom.get_vector()
            name = PDB_to_AC4DC_dict.get(atom.element)
            # Pop into our desired format
            if name == None or name not in allowed_atoms:
                pdb_atoms_ignored += atom.element + " "
                continue
            if name not in species_dict.keys():
                species_dict[name] = Atomic_Species(name,self) 
                pdb_atoms.append(atom.element)
            species_dict[name].add_atom(R)
        ac4dc_atoms_ignored = ""
        for string in allowed_atoms:
            if string not in species_dict.keys():
                ac4dc_atoms_ignored += string + " "
                continue
        print("The following atoms will be considered (AC4DC names ; pdb names):")
        for p,a in zip(pdb_atoms,[PDB_to_AC4DC_dict[x] for x in pdb_atoms]): 
            print("%-10s %10s" % (p, a))
        if pdb_atoms_ignored != "":
            print("The following pdb atoms were found but ignored:",pdb_atoms_ignored)
        if ac4dc_atoms_ignored != "":
            print("The following atoms were allowed but not found:",ac4dc_atoms_ignored)
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
    def __init__(self, experiment_name, photon_energy, detector_distance, hemisphere_screen = True, orientation_set = [[0,0,0],], x_orientations = 1, y_orientations = 1,q_cutoff="hemisphere", pixels_per_ring = 400, num_rings = 50,t_fineness=100):
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
            contains each set of cardan angles for each orientation to images. Crystal only. 
            Overriden for imagings with random orientations. TODO replace x_orientations and y_orientations with this.
        """
        self.experiment_name = experiment_name
        self.hemisphere_screen = hemisphere_screen
        self.photon_momentum = 2*np.pi/E_to_lamb(photon_energy)
        self.photon_energy = photon_energy
        self.detector_distance = detector_distance
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.x_orientations = x_orientations
        self.y_orientations = y_orientations

 
        if q_cutoff == "hemisphere":
            # as q = 2ksin(theta), non-inclusive upper bound is 45 degrees(theoretically).  
            q_cutoff = 2*self.photon_momentum/np.sqrt(2)  - 0.00000001       
        if q_cutoff == "sphere":
            q_cutoff = 2*self.photon_momentum  - 0.00000001
        self.q_cutoff = q_cutoff

        self.t_fineness = t_fineness

        self.phi_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)
        self.y_rotation = 0 # Current rotation of crystal (y axis currently)
        self.x_rotation = 0 # Current rotation of crystal (y axis currently)
        self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))        
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

        #TODO SPI needs to be fixed after being demolished by crystal-necessitated refactoring
        if SPI:
            result = Results_SPI()
            result.I = 0            
            q_samples = np.linspace(0,self.q_cutoff,self.num_rings)
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))      
                self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    for i, q in enumerate(q_samples):
                        ring[i] = self.generate_ring(q)
                        #ring[i].I = (ring[i].I+1)
                        #print("q:",q, "x:",ring[i].R,"I[alph=0]",ring[i].I[0])


                    # Initialise stuff that is constant between images (done here to access ring radii.)
                    if rot_y  == 0 and rot_x == 0:
                        phi = self.phi_array
                        #radii = np.zeros(self.num_rings)
                        #for i in range(len(ring)):
                            #radii[i] = ring[i].r           
                        #r, phi = np.meshgrid(radii, phi)     
                        q, phi = np.meshgrid(q_samples,phi)

                    I = np.zeros(q.shape)  # I is the intensity of the plot colour.
                    for ang in range(len(I)):
                        for pos in range(len(I[ang])):
                            I[ang][pos] = ring[pos].I[ang]                                         
                
                    result.I += I/(self.y_orientations*self.x_orientations)
                self.x_rotation += 2*np.pi/self.x_orientations
                result.phi = phi
                result.q = q
            return result

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
                #TODO vectorise
                for i, G in enumerate(bragg_points):   # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                    # if G[0] > largest_x:
                    #     largest_x = G[0]
                    #     print("New high G[0]!",G[0])
                    point = self.generate_point(G,cardan_angles)
                    #tmp_radii[i] = point.r
                    result.phi[i] = point.phi
                    result.phi_aligned[i] = point.phi_crystal_aligned
                    result.q[i] = point.q_parr_screen
                    result.I[i] = point.I
                    result.image_index[i] = j
                result.package_up(miller_indices)
                #Save the result object into its own file within the output folder for the experiment
                fpath = directory + str(cardan_angles) +".pickle"
                with open(fpath,"wb") as pickle_out:
                    pickle.dump(result,pickle_out)
            return used_orientations

    class Feature:
        def __init__(self,q,q_parr_screen,theta):
            self.q = q
            self.q_parr_screen = q_parr_screen # magnitude of q parallel to screen
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
        #r = self.q_to_r(q)# TODO
        theta = self.q_to_theta(q)
        q_parr_screen = self.q_to_q_scr(q)
        if self.hemisphere_screen:
            print("hemisphere screen not supported for SPI yet.")
        ring = self.Ring(q,q_parr_screen,theta)
        ring.I = self.illuminate(ring)
        return ring 
    def generate_point(self,G,cardan_angles): # G = vector
        q = np.sqrt(G[0]**2+G[1]**2+G[2]**2)
        #r = self.q_to_r(q) # TODO
        q_parr_screen = self.q_to_q_scr(q)
        if self.hemisphere_screen:
            q_parr_screen = self.q_to_q_scr_curved(G)
        theta = self.q_to_theta(q)
        point = self.Spot(q,q_parr_screen,theta)
        point.phi = np.arctan2(G[1],G[0])
        point.G = G
        #Crystal-aligned phi. i.e. we realign all the images to be in the 0 0 0 orientation.
        G = self.rotate_G_to_orientation(G,*cardan_angles,inverse=True)[0]
        point.phi_crystal_aligned = np.arctan2(G[1],G[0])
                
        if DEBUG:
            print("q_parr_screen",q_parr_screen,"phi",point.phi,"G",G[0],G[1],G[2])


        #print("phi",point.phi,"q_parr_screen",point.q_parr_screen)
        
        phis = np.array([point.phi])
        point.I = self.illuminate(point,phis,cardan_angles)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if q_parr_screen != check:
        #     print("Error, q_parr_screen =",q_parr_screen,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
    def illuminate(self,feature,phis = None,cardan_angles = None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet."""
        if phis == None:
            phis = self.phi_array
        F = np.zeros(phis.shape,dtype="complex_")
        count = 0
        for species in self.target.species_dict.values():
            species.set_scalar_form_factor(feature.q)
            # iterate through each symmetry of unit cell (each asymmetric unit)
            for s in range(len(self.target.sym_factors)):
                for R in species.coords:
                        count += 1
                        # Rotate to crystal's current orientation 
                        R = R.left_multiply(self.y_rot_matrix)  
                        R = R.left_multiply(self.x_rot_matrix)   
                        coord = np.multiply(R.get_array(),self.target.sym_factors[s]) + np.multiply(self.target.cell_dim,self.target.sym_translations[s])
                        # Get spatial factor T
                        T = np.zeros(phis.shape,dtype="complex_")
                        if SPI:
                            T = self.SPI_spatial_factor(phis,coord,feature)
                        else:
                            T= self.spatial_factor(phis,coord,feature,cardan_angles)
                        F += species.ff*T
                        # Rotate atom for next sample            
        #print("iterated through",count,"atoms")
        I = np.square(np.abs(F))
        

        return I 

    def spatial_factor(self,phi_array,coord,feature,cardan_angles):
        """ theta = scattering angle relative to z-y plane """ 
        q_z = feature.G[2]#feature.q_parr_screen*np.sin(feature.theta) TODO replace with functional method
        q_z = np.tile(q_z,len(phi_array))
        # phi angle of vector relative to x-y plane, pointing from screen centre to point hit. 
        q_y = feature.G[1]#q_z*np.sin(phi_array)
        q_x = feature.G[0]#q_z*np.cos(phi_array)
        q_vect = np.column_stack([q_x,q_y,q_z])
        # Rotate our q vector BACK to the 0 0 0 orientation of the crystal.
        q_vect = self.rotate_G_to_orientation(q_vect,*cardan_angles,inverse=True)[0]            
        spatial_structure_factor = np.exp(-1j*np.dot(q_vect,coord))   
        return spatial_structure_factor

    def SPI_spatial_factor(self,phi_array,R,feature):
        """ theta = scattering angle relative to z-y plane """ 
        q_z = feature.q*np.sin(feature.theta)
        q_z = np.tile(q_z,len(phi_array))
        # phi angle of vector relative to x-y plane, pointing from screen centre to point hit. 
        q_y = feature.q_parr_screen*np.sin(phi_array)
        q_x = feature.q_parr_screen*np.cos(phi_array)
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
        Currently only checked to work with cubics. 
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
        # TODO double check not cutting off possible values.
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
        indices = np.array([*set(indices)])   
    
        # Selection rules
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

        q_cutoff_rule = lambda g: np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2))<= self.q_cutoff
        #q_cutoff_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.q_cutoff
        q0 = self.photon_momentum
        
        # Catch the miller indices with a boolean mask
        mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(q_cutoff_rule,1,G_temp)*np.apply_along_axis(self.mosaic_elastic_condition,1,G_temp)
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

    def rotate_G_to_orientation(self,G,alpha=0,beta=0,gamma=0,random=False,inverse = False):
        '''
        Rotates the G vectors to correspond to the crystal when oriented to the given cardan angles. 
        alpha: x-axis angle [radians].
        beta:  y-axis angle [radians].
        gamma: z-axis angle [radians].
        ## Returns:
        G: transformed G
        list: [alpha,beta,gamma] (the angles used)
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
        rot_matrix = R.as_matrix()
        if inverse: 
            rot_matrix = rot_matrix.T # transpose == inverse
        # Apply rotation
        #print("Rotation matrix:")
        #print(R.as_matrix())
        G = np.matmul(np.around(rot_matrix,decimals=10),G)

        if transposed:
            G = G.T
        return G, [alpha,beta,gamma]
                
    # q_abs [1/bohr]
    def q_to_theta(self,q):
        lamb = E_to_lamb(self.photon_energy)
        theta = np.arcsin(lamb*q/(4*np.pi)) 
        #if not (theta >= 0 and theta < 45):
            #print("warning, theta out of bounds")        
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
            return q/(sin(np.pi/2-theta))#q*np.cos(theta)
    def q_to_q_scr_curved(self,G):
        return np.sqrt(G[0]**2+G[1]**2)
    
    def mosaic_elastic_condition(self,q_vect):
        ''' determines whether vector q is allowed.'''
        # Rocking angle version
        q = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + q_vect[2]**2)
        if q > self.q_cutoff:
            return False
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
        #
        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.

def scatter_scatter_plot(neutze_R = True, crystal_aligned_frame = False ,SPI_result = None,full_range = True,num_arcs = 50,num_subdivisions = 40, result_handle = None, compare_handle = None, normalise_intensity_map = False, show_grid = False, cmap_power = 1, cmap = None, min_alpha = 0.05, max_alpha = 1, bg_colour = "black",solid_colour = "white", show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False,  fixed_dot_size = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None):
    ''' (Complete spaghetti at this point.)
    Plots the simulated scattering image.
    result_handle:

    compare_handle:
        results/compare_handle/ is the directory of results that will be subtracted from those in the results directory.
        Must have same orientations.
    
    '''
    mpl.rcParams['text.color'] = "blue"
    mpl.rcParams['axes.labelcolor'] = "blue"
    mpl.rcParams['xtick.color'] = "black"
    mpl.rcParams['ytick.color'] = "blue"    
    #https://stackoverflow.com/questions/26108436/how-can-i-get-the-matplotlib-rgb-color-given-the-colormap-name-boundrynorm-an
    if cmap != None:
        class MplColorHelper:

            def __init__(self, cmap_name, start_val, stop_val):
                self.cmap_name = cmap_name
                self.cmap = plt.get_cmap(cmap_name)
                self.norm = mpl.colors.Normalize(vmin=start_val, vmax=stop_val)
                self.scalarMap = cm.ScalarMappable(norm=self.norm, cmap=self.cmap)

            def get_rgb(self, val):
                if compare_handle == None:
                    # Plot actual image!
                    val = ((val-min_z)/(max_z-min_z))**cmap_power
                if val < 0 or val > 1:
                    print("error in val!")
                return self.scalarMap.to_rgba(val)            
    else:
        # Broken TODO
        r,g,b = to_rgb(solid_colour)
        colours = [(r,g,b,a) for a in np.clip(z/max_z,min_alpha,max_alpha)]
    
    def add_screen_properties():
        if radial_lim:
            bottom,top = plt.ylim()
            plt.ylim(bottom,radial_lim)
        if log_radial:
            plt.yscale("log")  
        plt.gca().set_facecolor(bg_colour) 
        plt.gcf().set_figwidth(20)
        plt.gcf().set_figheight(20)             
    
    if result_handle != None:
        results_dir = "results/"+result_handle+"/"
        compare_dir = None
        if compare_handle!= None:
            compare_dir = "results/"+compare_handle+"/"        
        ## Point-like (Crystalline)
        # Initialise R factor sector comparison plot. 
        sector_histogram = np.zeros((num_arcs,num_subdivisions)).T
        sector_num_histogram = np.zeros(sector_histogram.shape)
        sector_den_histogram = np.zeros(sector_histogram.shape)
        R_histogram = np.zeros((num_arcs,num_subdivisions)).T
        R_num_histogram = np.zeros(R_histogram.shape)
        R_den_histogram = np.zeros(R_histogram.shape)
        # Bin edges
        phi_edges = np.linspace(-np.pi,np.pi,num_arcs+1)
        radial_edges = np.linspace(0,radial_lim,num_subdivisions+1)        
        sector_phi,sector_radial = np.meshgrid(phi_edges,radial_edges)
        num_orientations = 0
        def get_old_histogram_contribution(I,phi,radial_axis):
                # We don't set density = True, because we don't want to normalise the weights.
                num_samples = np.histogram2d(phi, radial_axis, bins=(phi_edges, radial_edges))[0]
                num_samples[num_samples == 0] = 1 # avoid division by zero
                H = np.histogram2d(phi, radial_axis, weights=I, bins=(phi_edges, radial_edges))[0]
                H = H.T
                H = np.divide(H, num_samples.T)
                return H/num_orientations   
        
        def get_sector_histogram_contribution(I_ideal,I_real,phi,radial_axis,phi_edges = phi_edges):
            '''I_ideal = Intensities of undamaged target
               I_real = Intensities of damaged target
            '''#                                                                            _|-|_
            # R = Σ|sqrt(I_ideal) - sqrt(I_real|) / (Σ sqrt(I_ideal))  |   me -> ( ._.)ヽ(￣┏＿┓￣ R)  "thousands of lines of code just for you Señor R factor".
            N = np.abs(np.sqrt(I_ideal) - np.sqrt(I_real))
            D = np.sqrt(I_ideal)
            numerators = np.histogram2d(phi, radial_axis, weights=N, bins=(phi_edges, radial_edges))[0]
            denominators  = np.histogram2d(phi, radial_axis, weights=D, bins=(phi_edges, radial_edges))[0]
            numerators = numerators.T; denominators = denominators.T
            return numerators,denominators


        def plot_sectors(sector_histogram):            
            plt.close()
            fig = plt.figure()
            ax2 = fig.add_subplot(projection="polar",aspect="equal")
            sector_histogram = np.ma.masked_where(sector_histogram ==0, sector_histogram)  # masking for colour (replace with black)
            if full_range:
                pcolour = ax2.pcolormesh(sector_phi,sector_radial,sector_histogram,cmap=cmap,vmin=0,vmax=1)
            else: pcolour = ax2.pcolormesh(sector_phi,sector_radial,sector_histogram,cmap=cmap)
            fig.colorbar(pcolour)
            # print(sector_phi)
            # print()
            # print(sector_radial)
            # print()
            # print(sector_histogram)
            add_screen_properties()       
            plt.show()        

        fig = plt.figure()
        ax = fig.add_subplot(projection="polar")

        if compare_handle != None:
            if log_I:
                print("Not plotting logarithmic I, not supported for comparisons.")
                log_I = False  
            
            all_I_ideal = []
            all_I_real = []
        
        # Iterate through each orientation (file) to get minimum/maximum for normalisation of plot of scattering image (not difference image).
        max_z = -np.inf
        min_z = np.inf
        for filename in os.listdir(results_dir):
            #result = get_result(filename)
            result1,result2 = get_result(filename,results_dir,compare_dir)
            if result1 == "__PASS__":
                continue            
            if result1 == None:
                break  
            num_orientations += 1
            if compare_dir != None:
                max_z = max(max_z,np.max(result1.R))
                min_z = min(min_z,np.min(result1.R))
                if max_z > 1:
                    print("error, max_R =",max_z)
                if min_z < 0:
                    print("error, min_R =",min_z)
            else:
                max_z = max(max_z,np.max(result1.I))
                min_z = min(min_z,np.min(result1.I))
        # Apply log (to non-diff image)
        if log_I:
            max_z = np.log(max_z)
            min_z = np.log(min_z)
        print("max,min",max_z,min_z)
        # Plot each orientation's scattering pattern/add each to sector histogram
        added_colorbar = False
        for filename in os.listdir(results_dir):
            result1,result2 = get_result(filename,results_dir,compare_dir)
            if result1 == "__PASS__":
                continue
            if result1 == None:
                break         
            # get points at same position on screen.
            radial_axis = result1.r
            if plot_against_q:
                radial_axis =  result1.q        
            radial_axis = radial_axis[0]    
            identical_count = np.zeros(result1.I.shape)
            
            
            if crystal_aligned_frame:
                phi = result1.phi_aligned
            else:
                phi  = result1.phi
            tmp_I1 = np.zeros(result1.I.shape)
            tmp_I2 = np.zeros(tmp_I1.shape)
            processed_copies = []
            unique_values_mask = np.zeros(result1.I.shape,dtype="bool")
            # Catch for multiple overlapping points (within same orientation only!!!) (does work, but would be unlikely.)
            # TODO need to do something about fact that overlapping points with many plots will hide points. Not critical rn thanks to sectors. 
            for i in range(len(result1.I)):
                if (radial_axis[i], phi[i],result1.image_index[i])  in processed_copies:
                    unique_values_mask[i] = False
                    continue
                else:
                    unique_values_mask[i] = True
                    
                # Boolean masks
                matching_rad_idx = np.nonzero(radial_axis == radial_axis[i])  # numpy note: equiv. to np.where(condition). Non-zero part irrelevant.
                matching_phi_idx = np.nonzero(phi == phi[i])
                for elem in matching_rad_idx[0]:
                    if elem in matching_phi_idx[0]:
                        identical_count[i] += 1
                        matching1 = result1.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                        for value in matching1:
                            tmp_I1[i] += value    
                        if result2 != None:
                            matching2 = result2.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                            for value in matching2:
                                tmp_I2[i] += value                                                      

                #processed_copies.append((radial_axis[i],phi[i],result.image_index[i]))
            #print("processed copies:",processed_copies)


            result = copy.deepcopy(result1)
            result.I = tmp_I1 
            # get z used for colour of scatter plot.
            if compare_dir != None:
                result2.I = tmp_I2                    
                result.diff(result2)
                z = result.R[unique_values_mask]
            else:
                z = result.I[unique_values_mask]

            identical_count = identical_count[unique_values_mask]
            # Intensities used for histogram 
            I1 = tmp_I1[unique_values_mask]
            I2 = tmp_I2[unique_values_mask]
            radial_axis = radial_axis[unique_values_mask]
            phi = phi[unique_values_mask]

            #sector_histogram += get_histogram_contribution(z,result.phi,radial_axis)
            if result2 != None:
                ## Sectors/Fragmented rings ##
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis)
                # # Average over orientations' sector R factors:  # Doesn't work because when sector is empty it reduces the average.
                # non_zero_denon = denominators.copy()
                # non_zero_denon[non_zero_denon == 0] = 1
                # sector_histogram += np.divide(numerators,non_zero_denon)/num_orientations
                # Alternative: Not averaging:
                sector_num_histogram += numerators
                sector_den_histogram += denominators                

                ## Full rings ## Note: we are filling 2d arrays with same-valued arcs for each radius so that pcolormesh can plot circles. 
                # The numerators/denominators are effectively tiled, e.g.: numerators = np.tile(numerators,(1,num_arcs))
                # Average over entire rings.
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis,phi_edges = np.array([-np.pi,np.pi]))
                # non_zero_denon = denominators.copy()
                # non_zero_denon[non_zero_denon == 0] = 1                
                # R_histogram += np.divide(numerators,non_zero_denon)/num_orientations # Doesn't work because when sector is empty it reduces the average.
                # Alternative: Not averaging:
                R_num_histogram += numerators
                R_den_histogram += denominators            

                if neutze_R:
                    all_I_ideal.extend(I1)
                    all_I_real.extend(I2)


            if log_I: 
                z = np.log(z)      

            #debug_mask = (0.01 < radial_axis[0])*(radial_axis[0] < 100)
            #print(identical_count[debug_mask])
            #print(radial_axis[0][debug_mask])
            #print(phi[debug_mask ]*180/np.pi)
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

            # use cmap to get colours but with alpha following a specific rule
            COL = MplColorHelper(cmap, 0, 1) 
            alpha_modified_cmap_colours = np.array([[COL.get_rgb(K)[0], COL.get_rgb(K)[1],COL.get_rgb(K)[2],np.clip((K*(max_alpha-min_alpha))/(max_z) + min_alpha,min_alpha,None)] for K in z])
            #print("THING",thing)
            colours = [(r,g,b,a) for r,g,b,a in alpha_modified_cmap_colours] 
            colours = np.around(colours,10)   

            if not added_colorbar:
                if compare_dir == None: 
                    print("Warning, colorbar not working with color power at present. Need to create cmap from COL")
                    if not normalise_intensity_map:
                        # Good for debugging
                        fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.min(z),vmax=np.max(z)),cmap=cmap),ax=ax)
                        print("Attention: Not normalising intensities, but still arbitrary units")
                        print("Warning: not yet taking into account combined dots!! Scale is off!")
                    else:
                        fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,vmax=1),cmap=cmap),ax=ax) 
                else:
                    #TODO get this truncation of bar working https://stackoverflow.com/questions/40982050/matplotlib-how-to-cut-the-unwanted-part-of-a-colorbar
                    # from matplotlib import colorbar
                    # colors = cmap(np.linspace(1.-(0.5-0.3)/float(0.5), 1, cmap.N))
                    # cbar_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(cmap, colors)
                    # cax,_ = colorbar.make_axes(ax)
                    # norm= mpl.colors.Normalize(vmin=0,vmax=1)
                    # cbar = colorbar.ColorbarBase(cax, cmap=cbar_cmap, norm=norm)
                    # cbar.set_ticks([0.3,0.4,0.5])
                    # cbar.set_ticklabels([0.3,0.4,0.5])
                    fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=0,vmax=1),cmap=cmap),ax=ax) 
                         
                #print(z)
                added_colorbar = True

            sc = ax.scatter(phi,radial_axis,c=colours,s=s) #TODO try using alpha = 0.7 or something to guarantee overlapped points not hidden
            plt.grid(alpha = min(show_grid,0.6),dashes=(5,10)) 
            #TODO only show first index or something. or have option to switch between image indexes. oof coding.
            if show_labels:
                for i, miller_indices in enumerate(result.miller_indices[unique_values_mask]):
                    ax.annotate("%.0f" % miller_indices[0]+","+"%.0f" % miller_indices[1]+","+"%.0f" % miller_indices[2]+",", (phi[i], radial_axis[i]),ha='center') 
        add_screen_properties()
        plt.show()  
        if compare_handle != None:
            #print("Plotting orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
            #plot_sectors(sector_histogram = sector_histogram)
            non_zero_denon_histogram = sector_den_histogram.copy()
            non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
            sector_histogram = np.divide(sector_num_histogram,non_zero_denon_histogram)
            print("Plotting total R factor")
            plot_sectors(sector_histogram = sector_histogram)

            #print("plotting full ring orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
            #plot_sectors(R_histogram)

            
            non_zero_denon_histogram = R_den_histogram.copy()
            non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
            R_histogram = np.divide(R_num_histogram,non_zero_denon_histogram)       
            print("plotting full ring total R factor ") 
            plot_sectors(R_histogram) 

            if neutze_R:
                print("Neutze R: ")
                sqrt_ideal = np.sqrt(all_I_ideal)
                sqrt_real = np.sqrt(all_I_real)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)
                #neutze_histogram = np.histogram2d(phi, radial_axis, weights=R, bins=(np.array([-np.pi,np.pi]), radial_edges))[0]             
                #plot_sectors(neutze_histogram)

        

       
       
       
       
       
    else:
        result = SPI_result
        ## Continuous (SPI)
        if log_I: 
            z = np.log(result.I)
        else:
            z = result.I        
        if log_I and cutoff_log_intensity != None:
            #cutoff_log_intensity = -1
            z -= cutoff_log_intensity
            z[z<0] = 0
            pass
        #radial_axis = result.r
        if plot_against_q:
            radial_axis =  result.q         
        fig = plt.figure()
        ax = fig.add_subplot(projection="polar")
        phi_mesh = result.phi_mesh
        if crystal_aligned_frame:
            phi_mesh = result.phi_aligned_mesh        
        ax.pcolormesh(phi_mesh, radial_axis, z,cmap=cmap)
        ax.plot(result.phi, radial_axis, color = 'k',ls='none')
        plt.grid(alpha=min(show_grid,0.6),dashes=(5,10)) 

        if radial_lim:
            bottom,top = plt.ylim()
            plt.ylim(bottom,radial_lim)
        if log_radial:
            plt.yscale("log")

    def get_orientation_set_of_folder():
        '''Returns a list of orientations present in results subfolder'''



# https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
from mpl_toolkits.axes_grid1 import AxesGrid
import matplotlib

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def get_result(filename,results_dir,compare_dir = None):
    #Requires all orientations of result_handle in compare_handle, but not vice versa.
    fpath = os.path.join(results_dir, filename)
    if os.path.isfile(fpath):
        with open(fpath,'rb') as f:
            result1 = pickle.load(f)
    else: 
        return "__PASS__" "__PASS__"
    result2 = None
    if compare_dir != None:
        if filename in os.listdir(compare_dir):
            fpath2 = os.path.join(compare_dir, filename)
            with open(fpath2,'rb') as f:
                result2 = pickle.load(f)     
                result1.diff(result2)             
        else:
            print("ERROR, missing matching orientation in comparison directory")
            return None, None # No corresponding file found.          
    return result1, result2




#%% https://scripts.iucr.org/cgi-bin/paper?S0021889807029238, http://superflip.fzu.cz/

import pandas as pd
import csv
def create_reflection_file(result_handle,overwrite=False):
    print("Creating reflection file for",result_handle)
    directory = "reflections/"    
    results_dir = "results/"+ result_handle+"/"
    os.makedirs(directory, exist_ok=True) 
    results_dir = "results/"+result_handle+"/"
    init = False
    for filename in os.listdir(results_dir):
        result = get_result(filename,results_dir)[0]
        if result == "__PASS__":
            continue            
        if result == None:
            break          
        miller_indices = result.miller_indices
        intensity = np.array([result.I]).T
        data = np.concatenate((miller_indices,intensity),axis=1)
        # Put in data frame
        columns = ["h","k","l","I"]
        new_df = pd.DataFrame(data=data, columns=columns)
        if init == False:
            df = new_df
            init = True
        else:
            df = pd.concat((df,new_df))
    # Save as file
    out_path = directory + result_handle + ".rfl"
    if os.path.isfile(out_path): 
        if not overwrite:
            print("Cannot write, file already present at",out_path)
            return
        os.remove(out_path)
    columns.reverse()
    df = df.sort_values(by=["l","k","h"],axis=0)
    for i in ("hkl"):
        df[i] = df[i].astype('int')
    df["I"] = df["I"].astype('float')
    df = df.round(6)
    df.drop_duplicates(inplace=True) 
    df.to_csv(out_path,header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")

#create_reflection_file("f1_11",True)

#%% 
root = "lys"
tag = 0
#%%
### Simulate tetrapeptide
tag += 1
#  #TODO make this nicer and pop in a function 
DEBUG = False
energy = 12561 # Tetrapeptide (though AC4DCsimulation was 6000 eV oops) 
#energy =17445   #crambin - from resolutions. in pdb file, need to double check calcs: #q_min=0.11,q_maxres = 3.9,  my defaults: pixels_per_ring = 400, num_rings = 200


exp_name1 = str(root) + "1_" + str(tag)
exp_name2 = str(root) + "2_" + str(tag)



allowed_atoms_1 = ["C_fast","N_fast","O_fast","Sodion"]#,"S_fast"]
allowed_atoms_2 = ["C_fast","N_fast","O_fast"]#,"S_fast"]

end_time_1 = -10#-9.95
end_time_2 = -10#0#-9.80  

#orientation_set = [[5.532278012665244, 1.6378991611963682, 1.3552062099726534]] # Spooky spider orientation
num_orients = 100#50
# [ax_x,ax_y,ax_z] = vector parallel to axis. Overridden if random orientations.
ax_x = 0
ax_y = 0
ax_z = 0
random_orientation = True # if true, overrides orientation_set

rock_angle = 0.3 # degrees

pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb" #tetrapeptide

output_handle = "sodion_tetrapeptide_2"

hemisphere_screen = True

#---------------------------------#


# if not random:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

experiment1 = XFEL(exp_name1,energy,100, hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)
experiment2 = XFEL(exp_name2,energy,100, hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)


#0.85  - 0.88 angstrom resolution - >  

# sym_translations = [np.array([0,0,0])]
# cell_dim = [np.array([1,1,1])]  

crystal = Crystal(pdb_path,allowed_atoms_1,rocking_angle = rock_angle*np.pi/180,CNO_to_N=False,cell_packing = "SC")
crystal.set_cell_dim(4.813 ,  17.151 ,  29.564)
crystal.add_symmetry(np.array([-1, -1,1]),np.array([0.5,0,0.5]))  #2555
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0.5]))  #3555
crystal.add_symmetry(np.array([1, -1,-1]),np.array([0.5,0.5,0]))  #4555

SPI = False

exp1_orientations = experiment1.spooky_laser(-10,end_time_1,output_handle,crystal, random_orientation = random_orientation, SPI=SPI)
create_reflection_file(exp_name1)
experiment2.orientation_set = exp1_orientations
experiment2.spooky_laser(-10,end_time_2,output_handle,crystal, random_orientation = False, SPI=SPI)

stylin()
# %%
# Plot atoms retrieved from Superflip/EDMA
df = pd.read_csv('tet_plot_points.pl',delim_whitespace=True)
x = df['x']
y = df['y']
z = df['z']


xscale = 4.813 
yscale = 17.151
zscale = 29.564
fig =plt.figure(figsize =(zscale,yscale))
#ax = Axes3D(fig)
#ax.scatter(x,y,z)


x*= xscale
y*= yscale
plt.scatter(z,y,c=-x,cmap="Blues",s=200)
plt.xlabel("z (Ang)")
plt.ylabel("y (Ang)")

#%%
#%%
### Simulate lysozyme
root = "lys"
tag = "v1"
#  #TODO make this nicer and pop in a function 
DEBUG = False
energy = 6000 #


exp_name1 = str(root) + "1_" + str(tag)
exp_name2 = str(root) + "2_" + str(tag)

CNO_to_N = True

allowed_atoms_1 = ["N_fast","S_fast"]
allowed_atoms_2 = ["N_fast","S_fast"]

end_time_1 = -10#-9.95
end_time_2 = 10#0#-9.80  

num_orients = 50
# [ax_x,ax_y,ax_z] = vector parallel to axis. Overridden if random orientations.
ax_x = 1
ax_y = 1
ax_z = 0
random_orientation = True # if true, overrides orientation_set

rock_angle = 0.3 # degrees

pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb" #tetrapeptide

output_handle = "D_lys_neutze_simple_7"

hemisphere_screen = True

#---------------------------------#


# if not random:
if ax_x == ax_y == ax_z and ax_z == 0 and random_orientation == False:
    throw
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

experiment1 = XFEL(exp_name1,energy,100, hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)
experiment2 = XFEL(exp_name2,energy,100, hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)


#0.85  - 0.88 angstrom resolution - >  

# sym_translations = [np.array([0,0,0])]
# cell_dim = [np.array([1,1,1])]  

crystal = Crystal(pdb_path,allowed_atoms_1,rocking_angle = rock_angle*np.pi/180,CNO_to_N=CNO_to_N,cell_packing = "SC")
crystal.set_cell_dim(79.000  , 79.000  , 38.000)
crystal.add_symmetry(np.array([-1, -1,1]),np.array([0.5,0,0.5]))  #2555
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0.5]))  #3555
crystal.add_symmetry(np.array([1, -1,-1]),np.array([0.5,0.5,0]))  #4555

SPI = False

exp1_orientations = experiment1.spooky_laser(-10,end_time_1,output_handle,crystal, random_orientation = random_orientation, SPI=SPI)
create_reflection_file(exp_name1)
experiment2.orientation_set = exp1_orientations
experiment2.spooky_laser(-10,end_time_2,output_handle,crystal, random_orientation = False, SPI=SPI)

stylin()

#%%
### Simulate crambin
tag += 1
#  #TODO make this nicer and pop in a function 
DEBUG = False
energy =17445   #crambin - from resolutions. in pdb file, need to double check calcs: #q_min=0.11,q_cmaxres = 3.9,  my defaults: pixels_per_ring = 400, num_rings = 200


exp_name1 = str(root) + "1_" + str(tag)
exp_name2 = str(root) + "2_" + str(tag)



allowed_atoms_1 = ["C_fast","N_fast","O_fast","S_fast"]
allowed_atoms_2 = ["C_fast","N_fast","O_fast","S_fast"]

end_time_1 = -10#-9.95
end_time_2 = -10#-9.80  

#orientation_set = [[5.532278012665244, 1.6378991611963682, 1.3552062099726534]] # Spooky spider orientation
num_orients = 30
# [ax_x,ax_y,ax_z] = vector parallel to axis. Overridden if random orientations.
ax_x = 0
ax_y = 0
ax_z = 1
random_orientation = True  # if true, overrides orientation_set

rock_angle = 0.15 # degrees
hemisphere_screen = True

pdb_path = "/home/speno/AC4DC/scripts/scattering/3u7t.pdb" #Crambin

output_handle = "Improved_Lys_mid_6"

q_cutoff = "hemisphere" #"max"

test_the_onion = False  # <- a test case I'm familiar with, I'm using it to check if anything breaks.
#---------------------------------#


# if not random:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

if test_the_onion:
    orientation_set = [[0,0,0],[np.pi/4,0,0]]
    #q_cutoff = 10

experiment1 = XFEL(exp_name1,energy,100, q_cutoff = q_cutoff,hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)
experiment2 = XFEL(exp_name2,energy,100, q_cutoff = q_cutoff,hemisphere_screen = hemisphere_screen, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)


crystal = Crystal(pdb_path,allowed_atoms_1,rocking_angle = rock_angle*np.pi/180,CNO_to_N=False,cell_packing = "SC")
crystal.set_cell_dim(22.795 ,  18.826 ,  41.042)
if test_the_onion:
   crystal.set_cell_dim(20,20,20)
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0]))   #2555


SPI = False

exp1_orientations = experiment1.spooky_laser(-10,end_time_1,output_handle,crystal, random_orientation = random_orientation, SPI=SPI)
create_reflection_file(exp_name1)
experiment2.orientation_set = exp1_orientations
experiment2.spooky_laser(-10,end_time_2,output_handle,crystal, random_orientation = False, SPI=SPI)

stylin()
#

 
#%%
# stylin' 
def stylin():
    experiment1_name = exp_name1#"Lys_9.95_random"#exp_name1
    experiment2_name = exp_name2#"lys_9.80_random"#exp_name2 

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
    cmap = shiftedColorMap(matplotlib.cm.RdYlGn_r,midpoint=0.2)#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
    cmap.set_bad(color='black')
    cmap_power = 1.6
    min_alpha = 0.3
    max_alpha = 1
    colour = "y"
    radial_lim = 5
    full_crange_sectors = False

    cmap_intensity = "inferno"


    # screen_radius = 150#55#165    #
    # q_scr_lim = experiment.r_to_q(screen_radius) #experiment.r_to_q_scr(screen_radius)#3.9  #NOTE this won't be the actual max q_parr_screen but ah well.
    # zoom_to_fit = True
    # ####### n'
    # # plottin'

    # as q = ksin(theta).  
    q_scr_max = experiment1.q_cutoff#experiment1.q_cutoff*(1/np.sqrt(2)) # for flat screen. as q_scr = qcos(theta). (q_z = qsin(theta) =ksin^2(theta)), max theta is 45. (Though experimentally ~ 22 as of HR paper)

    zoom_to_fit = False
    if not zoom_to_fit:
        if use_q:
            radial_lim = q_scr_max+0.2
        #els:
    #else:
    #     radial_lim = screen_radius#min(screen_radius,experiment.q_to_r(experiment.q_cutoff))
    #     print(radial_lim)
    #     if use_q:
    #         #radial_lim = experiment.r_to_q_scr(radial_lim)
    #         radial_lim = q_scr_lim #radial_lim = min(q_scr_lim,experiment.q_to_q_scr(experiment.q_cutoff))
    #         print(radial_lim)
    # else:
    #     radial_lim = None

    #TODO fix above to work with distance

    # Sectors
    scatter_scatter_plot(crystal_aligned_frame = False,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
    # Intensity of experiment 1. 
    scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
    # Sectors
    scatter_scatter_plot(crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
    # Intensity of experiment 1. 
    scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)

    #NEED TO CHECK. We have a 1:1 mapping from q to q_parr, but with our miller indices we are generating multiple q_parr with diff q.
    # So we SHOULD get the same q_parr with different q_z. Which makes sense since we are just doing cosine. But still icky maybe?
    # Need to double check we get different intensities for same q_parr. Pretty sure that's implemented.
stylin()

#TODO 
# log doesnt work atm.
# Get the rings to correspond to actual rings
# %%
#%%
### Get undamaged lysozyme SPI
tag += 1
#  #TODO make this nicer and pop in a function 
DEBUG = False
energy = 12000 

allowed_atoms= ["C_fast","N_fast","O_fast","S_fast"]

end_time = -10#-9.95

rock_angle = 0.3 # degrees

pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb" #Lysozyme

output_handle = "Improved_Lys_mid_6"

#---------------------------------#


# if not random:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

experiment = XFEL("SPI",energy,100, hemisphere_screen = True, pixels_per_ring = 500, num_rings = 500,t_fineness=100)
crystal = Crystal(pdb_path,allowed_atoms,CNO_to_N=False,cell_packing = "SC")

SPI = True

SPI_result = experiment.spooky_laser(-10,end_time,output_handle,crystal, SPI=SPI)

#%%
#TODO Get the colors to match neutze.
use_q = True
log_radial = False
log_I = True
cutoff_log_intensity = None 
cmap = "PiYG_r"# "plasma"
radial_lim = 2
scatter_scatter_plot(SPI_result,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)

# get sinc functions with non-infinite crystal
# 
# %%
import numpy as np
theta = np.arcsin(10/(4*np.pi))
#Assuming theta > 0:
delt = 1
min_err = 5*np.sin(theta - delt/2)
max_err = 5*np.sin(theta+delt/2)
err = np.sqrt(5)**2
min_err <= err <= max_err
# %%
