#%%
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

#TODO make class

####
%colors nocolor
import os
os.getcwd()
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
print(sys.path)
######

from Bio.PDB.vectors import Vector as Bio_Vect
from Bio.PDB.vectors import homog_trans_mtx, set_homog_trans_mtx
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
import colorcet as cc; import cmasher as cmr
from mpl_toolkits.mplot3d import Axes3D

%matplotlib widget
plt.ioff()  # stops weird vscode stuff

DEBUG = False 
c_au = 137.036; eV_per_Ha = 27.211385; ang_per_bohr = 1/1.88973  # > 1 ang = 1.88973 bohr

class Results():
    def __init__(self,num_points,image_index):
        self.phi = np.zeros(num_points)
        self.phi_aligned = np.zeros(num_points)
        self.I = np.zeros(num_points)
        self.q = np.zeros(num_points)
        self.r = np.zeros(num_points)
        self.image_index = image_index 

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
class Results_Grid():
    pass    

class Crystal():
    def __init__(self, pdb_fpath, allowed_atoms, cell_dim, positional_stdv = 0, is_damaged=True, include_symmetries = True, rocking_angle = 0.3, cell_packing = "SC", CNO_to_N = False, orbitals_as_shells = True,cell_scale = 1):
        '''
        rocking_angle [degrees]
        cell_packing ("SC","BCC","FCC","FCC-D")
        cell_dim [angstroms]
        '''
        self.cell_packing = cell_packing
        self.rocking_angle = rocking_angle * np.pi/180
        self.cell_dim = cell_dim/ang_per_bohr  #np.array([1,1,1],dtype="float")              
        self.is_damaged = is_damaged
        self.cell_scale = cell_scale  # TODO allow for non-cubic crystals and non SC cell packing.
        self.pdb_fpath = pdb_fpath
        self.positional_stdv = positional_stdv/ang_per_bohr # RMS error in coord positions, SPI sim only.
        
        # Symmetry for each asymmetric unit simulated. If non-SPI, this is defining the supercell.
        self.sym_rotations = []; self.sym_translations = []; 
        if include_symmetries:
            self.parse_symmetry_from_pdb() # All asymmetric units in unit cell
        else:   
            self.add_symmetry(np.identity(3),np.zeros((3)))  # Only one asymmetric unit per unit cell

        parser=PDBParser(PERMISSIVE=1)
        structure_id = os.path.basename(self.pdb_fpath)
        structure = parser.get_structure(structure_id, self.pdb_fpath)
        
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

        for species in species_dict.values():
            species.set_coord_deviation()
        self.species_dict = species_dict 


    def set_ff_calculator(self,ff_calculator):
        self.ff_calculator = ff_calculator                  
    
    def add_symmetry(self,symmetry_factor,symmetry_translation,symmetry_label="Symmetry Operation:"):
        symmetry_translation/= ang_per_bohr
        # Simple cubic packing.
        if self.cell_packing == "SC":
            self.num_cells = self.cell_scale**3
            # Generate the coordinates of the cube (performance: defining scaling matrix/cube coords outside of here would be more efficient). But only runs once  \_( '_')_/ ¯\_(ツ)_/¯.
            x, y, z= np.meshgrid(np.arange(0, self.cell_scale), np.arange(0, self.cell_scale), np.arange(0, self.cell_scale))
            cube_coords = np.stack([ y.flatten(), x.flatten(), z.flatten()], axis = -1)
            # Construct the cube by adding the "unit cell translation" to the symmetry translation. 
            # (e.g. if crystal is centred on origin, in fractional crystallographic coordinates the index of the unit cell is equal to the unit cell's translation from the origin.) 
            print_thing = np.array(["][x]   [","][y] + [","][z]   ["],dtype=object)
            print(symmetry_label,"\n",np.c_[symmetry_factor, print_thing ,symmetry_translation],sep="")
            for coord in cube_coords:
                #print(coord)
                translation = coord*self.cell_dim + symmetry_translation 
                self.sym_translations.append(translation)
                self.sym_rotations.append(symmetry_factor)
                 
        else:
            #Can just duplicate symmetries if not SPI since equivalent for crystal.
            raise Exception("Lacking implementation")
            
    def parse_symmetry_from_pdb(target):
        '''
        Parses the symmetry transformations from the pdb file into ndarrays and adds each to the crystal.  
        Also returns the parsed matrices. 
        '''
        at_SYMOP = False
        at_symmetry_xformations = False
        sym_factor = np.zeros((3,3))
        sym_trans = np.zeros((3))
        sym_labels = []
        sym_mtces_parsed = []        
        with open(target.pdb_fpath) as pdb_file:
            for line in pdb_file:
                line = line.strip()
                # End data - Check if left section of symmetry data, marked by the line containing solely "REMARK 290". 
                if len(line) == 10: 
                    if at_SYMOP:
                        at_SYMOP = False
                    if at_symmetry_xformations: 
                        # Reached end of data
                        break
                # Symmetry operators (purely for terminal output - doesn't affect simulation)
                if at_SYMOP:
                    sym_labels.append(line[17:21]+": ("+line[21:].strip()+")")
                # Add each symmetry to target
                if at_symmetry_xformations:
                    entries = line.split()[2:]
                    x = int(entries[0][-1]) - 1  
                    entries = [float(f) for f in entries[2:]]
                    # x/y index rows/cols
                    for y, element in enumerate(entries[:-1]):
                        sym_factor[x,y] = element
                    sym_trans[x] = entries[-1]
                    if x == 2:
                        target.add_symmetry(sym_factor.copy(),sym_trans.copy(),sym_labels.pop(0))
                        sym_mtces_parsed.append([sym_factor.copy(),sym_trans.copy()])
                        
                
                # Start data - Check if this line marks the next as the beginning of the desired data.
                if "NNNMMM   OPERATOR" in line:
                    at_SYMOP = True
                    print("Parsing symmetry operations...")
                if line == "REMARK 290 RELATED MOLECULES.":
                    at_symmetry_xformations = True
        return sym_mtces_parsed
    
    def get_sym_xfmed_point(self,R,symmetry_index):
        '''
        R = 1D or 2D array of shape (3,N)
        '''
        i = symmetry_index
        # Transpose because we've stored the vectors in the 0th axis. 
        return (self.sym_rotations[i] @ R.T).T+ self.sym_translations[i]   # dim = [xyz,xyz] x [xyz,N] or dim = [xyz,xyz] x [xyz,1]
            
    def plot_me(self,max_points = 500):
        '''
        plots the symmetry points. Origin is the centre of the base asymmetric unit (not the centre of a unit cell).
        RMS error in positions ('positional_stdv') is not accounted for
        '''
        plt.close()
        fig = plt.figure(figsize=(4,4))
        ax = fig.add_subplot(111, projection='3d')
        num_atoms_avail = 0
        for species in self.species_dict.values():
             num_atoms_avail += len(species.coords)
        num_test_points = max(1,int(max_points/(self.num_cells*len(self.sym_translations))))
        num_test_points = min(num_atoms_avail,num_test_points)         
        test_points = np.empty((num_test_points,3))
        qualifier = "sample of"
        if num_atoms_avail == len(test_points):
            qualifier = "all" 
        print("Generating plot with",qualifier,len(test_points),"atoms per asymmetric unit")
        i = 0
        for species in self.species_dict.values():
            if i >= num_test_points:
                break                
            for R in species.coords:
                test_points[i] = R
                i+=1
                if i >= num_test_points:
                    break      
        plot_coords = []
        for i in range(len(self.sym_rotations)):
            coord_list = self.get_sym_xfmed_point(test_points,i).tolist()
            plot_coords.extend(coord_list)  
        color = np.empty(len(plot_coords),dtype=object) 
        color.fill('b')

        for i in range(int(len(self.sym_rotations)/self.num_cells)):           
            color[i*self.num_cells*len(test_points):i*self.num_cells*len(test_points)+len(test_points)] = 'c'  # atoms in one unit cell  
        color[:len(test_points)] = 'g'  # atoms in one asym unit       
        for i in range(self.num_cells):
            for j in range(len(self.sym_rotations)):
                color[i*self.num_cells + len(test_points)*j:i*self.num_cells + len(test_points)*j+ 1] = 'y'      # same atom in every asym unit          
            color[i*len(test_points):i*len(test_points) + 1] = 'r'      # same atom in every same unit cell         
        
        plot_coords = np.array(plot_coords)*ang_per_bohr
        # if np.unique(plot_coords) != plot_coords:
        #     a, unique_indices = np.apply_along_axis(np.unique,1,plot_coords,return_index=True)
        #     print(len(np.delete(plot_coords, unique_indices)), "non-unique coords found:")
        #     print(np.delete(plot_coords, unique_indices))
        coord_errors = 0  
        for i, coord in enumerate(plot_coords):
            for j, coord2 in enumerate(plot_coords):
                if coord[0] == coord2[0] and coord[1] == coord2[1]  and coord[2] == coord2[2]  and i!=j:
                    if coord_errors < 50:
                        print("error duplicate coord found:",coord)
                    
                    coord_errors += 1
        if coord_errors >=50:
            print("----------")
            print(coord_errors - 49,"more duplicate coord errors suppressed")
            print("----------")
                        
        #color[color!='y'] = '#0f0f0f00'   # example to see just one atom type
        print("---")
        print("x_min/max:",np.min(plot_coords[:,0]),np.max(plot_coords[:,0]))
        print("y_min/max:",np.min(plot_coords[:,1]),np.max(plot_coords[:,1]))
        print("z_min/max:",np.min(plot_coords[:,2]),np.max(plot_coords[:,2]))
        print("---")
        ax.scatter(plot_coords[:,0],plot_coords[:,1],plot_coords[:,2],c=color)
        ax.set_xlabel('x [$\AA$]')
        ax.set_ylabel('y [$\AA$]')
        ax.set_zlabel('z [$\AA$]')
        plt.show()

class Atomic_Species():
    def __init__(self,name,crystal):
        self.name = name 
        self.crystal = crystal 
        self.ff = 0 # form_factor
        self.coords = []  # coord of each atom in species in asymmetric unit

    def add_atom(self,vector):
        '''
        This function adds an atom to the asymmetric unit of the crystal. 
        We do not store additional coordinates, instead storing the symmetries, and an array of atomic states corresponding to each atom, for each symmetry. (so num symmetries * num atoms added)
        '''           

        self.coords.append(vector.get_array()/ang_per_bohr)
        
    def set_stochastic_states(self):
        '''
        We set a state for each atom, including symmetries, so that different q applied to the same atom at the same time corresponds to the same state. 
        When we are finding the atomic form factors, we call get_stochastic_f() on these states.
        The dimensionless nature of the model forces us to make the dubious approximation that an atom's state is independent of its prior states.
        With a hybrid molecular dynamics model informed by AC4DC, the nuclei's states could be tracked properly throughout time, and this function would be replaced
        by a call to the data of the atomic nuclei's states.
        '''
        self.times_used = self.crystal.ff_calculator.f_undamaged(0,self.name)[1]
        if self.num_atoms != len(self.crystal.sym_rotations)*len(self.coords):
            raise Exception("num atoms was not same on set_stochastic_states call as when set by set_coord_deviation")
        if self.crystal.is_damaged:
            # Initialise an array that tracks the individual atoms' form factors.
            orb_occs_shape = (self.num_atoms,) + self.crystal.ff_calculator.random_state_snapshots(self.name)[0].shape
            self.orb_occs= np.empty(orb_occs_shape)    # self.orb_occs[i] is an array of states corresponding to each time. We make the necessary approximation that an atom's state is independent of its prior states. This approximation is dubious at low unit cell numbers, but at higher numbers, because the contribution from an atom at the same relative cell coordinate and state as another atom will be equivalent, we get the same outcome so long as the probability distribution of states is representative of the actual distribution of states. i.e. tracking state history at the same global position is redundant at high unit cell counts where we can expect the distribution of states at a given coordinate to have a low deviation between species.  
            for idx in range(self.num_atoms):
                self.orb_occs[idx] = self.crystal.ff_calculator.random_state_snapshots(self.name)[0] 
    def set_coord_deviation(self):
        self.num_atoms = len(self.crystal.sym_rotations)*len(self.coords)
        self.error = np.empty((self.num_atoms,3))
        for idx in range(self.num_atoms):
            # get random error in spherical coords based on RMS error in position, convert to cartesian.
            # there's probably a better way to do it
            err_phi,err_thet = np.random.random()*2*np.pi, np.random.random()*np.pi
            err_r = np.random.normal(0,self.crystal.positional_stdv, size = (3))
            self.error[idx] = err_r*(np.sin(err_phi)*np.cos(err_thet),np.sin(err_phi)*np.cos(err_thet),np.cos(err_thet))


    def get_stochastic_f(atom_idx,q_arr):
        '''
        (Is defined by set_scalar_form_factor).
        Returns the form factor multiplied by sqrt(I). f.shape = ( len(times) , ) + momenta.shape  
        
        atom_idx, int or int array
        q_arr = mom. transfer [1/a0], scalar or array
        '''
        raise Exception("Did not set_scalar_form_factor before calling stochastic f")
        
    def set_scalar_form_factor(self,stochastic=True):
        '''
        
        '''
        # F_i(q) = f(q)*T_i(q), where f = self.ff is the time-integrated average 
        if not stochastic:
            pass
            #self.ff = self.crystal.ff_calculator.f_average(q_arr,self.name)      # note that in taking the integral to get this ff, we included the relative intensity.
        else:
            # Undamaged case, no stochastic dynamics.
            if not self.crystal.is_damaged: 
                def stoch_f_placeholder(atom_idx,q_arr): 
                    return self.crystal.ff_calculator.f_undamaged(q_arr,self.name)[0]
            # Damaged, we 
            else:
                def stoch_f_placeholder(atom_idx,q_arr): 
                    return self.crystal.ff_calculator.random_states_to_f_snapshots(self.times_used,self.orb_occs[atom_idx],q_arr,self.name)[0]  # f has form  [times,momenta]
            self.get_stochastic_f = stoch_f_placeholder
class XFEL():
    def __init__(self, experiment_name, photon_energy, detector_distance_mm=100, q_minimum = None, q_cutoff = None, max_triple_miller_idx = None, screen_type = "hemisphere", orientation_set = [[0,0,0],], x_orientations = 1, y_orientations = 1, pixels_per_ring = 400, num_rings = 50,t_fineness=100,SPI_y_rotation = 0,SPI_x_rotation = 0,SPI_z_rotation = 0):
        """ #### Initialise the imaging experiment's controlled parameters
        experiment_name:
            String that the output folder will be named.        
        photon_energy [eV]:
            Should be the same as that given in the original input file!
        detector_distance_mm [mm];
            The distance in mm between the target and the centre of the detector 
        screen_type:
            "circle", "hemisphere", "sphere". Circle corresponds to a flat screen - i.e. it 'squashes ya dots'. Sphere makes no difference if a max q is specified below the hemisphere q range.
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
        max_triple_miller_idx:
            If specified, the maximum momentum transfer is set to that corresponding to (m,m,m), where m = max_triple_miller_idx
        q_cutoff [1/angstrom]:
            If specified, only bragg points corresponding to momentum transfers at or below this value will be simulated.    
        """
        self.experiment_name = experiment_name
        self.detector_distance = detector_distance_mm*10e7/ang_per_bohr  # [a0 (bohr)]
        self.photon_momentum = 2*np.pi/E_to_lamb(photon_energy)  # atomic units, [a0^-1]. 
        self.photon_energy = photon_energy  #eV Attention: not in atomic units
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.x_orientations = x_orientations
        self.y_orientations = y_orientations
        self.max_triple_miller_idx = max_triple_miller_idx


        self.min_q = 0
        if q_minimum!=None:
            self.min_q = q_minimum/1.88973

        self.hemisphere_screen = True
        if screen_type == "flat": #well a circle
            self.hemisphere_screen = False
            self.max_q = self.photon_momentum - 0.00000001
        elif screen_type == "hemisphere":
            # as q = 2ksin(theta), non-inclusive upper bound is 45 degrees(theoretically).  
            self.max_q = 2*self.photon_momentum/np.sqrt(2)  - 0.00000001

        elif screen_type == "sphere":
            self.max_q = 2*self.photon_momentum  - 0.00000001
        else:
            raise Exception("Available screen types are 'flat, 'hemisphere', and 'sphere'")
        if q_cutoff != None:
            self.max_q = min(self.max_q,q_cutoff/1.88973) # convert to atomic units 


        self.t_fineness = t_fineness

        self.phi_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)  # used only for SPI. Really need to refactor...
        self.z_rotation = SPI_z_rotation *np.pi/180
        self.y_rotation = SPI_y_rotation *np.pi/180 # Current rotation of crystal (y axis currently) (did I mean z axis?)
        self.x_rotation = SPI_x_rotation *np.pi/180# Current rotation of crystal (y axis currently)
        self.z_rot_matrix = rotaxis2m(self.z_rotation,Bio_Vect(0, 0, 1))
        self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))        
        self.orientation_set = orientation_set 
                      
    def get_ff_calculator(self,start_time,end_time,damage_output_handle):
        ff_calculator = Plotter(damage_output_handle,"y")
        plt.close()
        ff_calculator.initialise_coherence_params(start_time,end_time,self.max_q,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return ff_calculator
    
    def spooky_laser(self, start_time, end_time, damage_output_handle, target, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False):
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
        for species in target.species_dict.values():
            species.set_stochastic_states()      
        self.target = target

        

        if pixels_across == None and circle_grid == False and SPI:
            raise Exception("Require pixels_across argument for rectangular screen")
        
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

        
        if SPI and circle_grid:
            ring = np.empty(self.num_rings,dtype="object")
            result = Results_SPI()
            result.I = 0            
            q_sep = (self.max_q-self.min_q)/(self.num_rings)
            q_samples = np.linspace(self.min_q,self.max_q,self.num_rings) + q_sep/2
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))      
                self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    for i, q in enumerate(q_samples):
                        # approximate angle subtended
                        angle_subtended = self.q_to_theta(q+q_sep/2) - self.q_to_theta(q-q_sep/2)
                        ring[i] = self.generate_ring(q,self.phi_array,angle_subtended)
                        #print("q:",q, "x:",ring[i].R,"I[alph=0]",ring[i].I[0])


                    # Initialise stuff that is constant between images (done here to access ring radii.)
                    if rot_y  == 0 and rot_x == 0:
                        phi = self.phi_array
                        #radii = np.zeros(self.num_rings)
                        #for i in range(len(ring)):
                            #radii[i] = ring[i].r           
                        #r, phi = np.meshgrid(radii, phi)     
                        q, phi = np.meshgrid(q_samples,phi)

                    I = np.zeros(q.shape)
                    for ang in range(len(I)):
                        for pos in range(len(I[ang])):
                            I[ang][pos] = ring[pos].I[ang]                                         
                
                    result.I += I/(self.y_orientations*self.x_orientations)
                self.x_rotation += 2*np.pi/self.x_orientations
                result.phi = phi
                result.q = q
            #result.package_up()
            return result
        elif SPI:
            def resolution_to_q_scr(d):
                q = 2*np.pi/d
                return self.q_to_q_scr(q)
            N = pixels_across#100  NxN cells
            cell = np.zeros((N,N),dtype="object")
            result = Results_Grid()      
            
            ### Geometry
            # In crystallography, for a resolution d we have q = 2*pi/d as lambda=2dsin(theta), q = (4pi/lambda)sin(theta). Possibly a questionable definition for SPI without periodicity, but it is used for consistency, and Nuetze 2000 makes no distinction.
            d = 2/ang_per_bohr # resolution            
            max_theta = self.q_to_theta(2*np.pi/d)
            #screen_width = resolution_to_X(d) * 2
            screen_distance = self.detector_distance # screen-target separation [a0]
            # res. at edge corner or centre? Surely at centre, for a full ring of information
            screen_width = 2*np.sin(2*max_theta)*screen_distance  # edge centre  # We have placed our screen to have the desired resolution at the "rim".
            #screen_width=  2 * (np.sin(2*max_theta)/np.sqrt(2)) * screen_distance # corner
            print("Screen width:",round(screen_width*ang_per_bohr/10e7,2),"mm")
            # X is the distance from the centre of the screen to the point of incidence (flat screen)
            # We know where the cells are, they are equally spaced. We also know the distance from the screen to the detector.
            # This gives us theta.
            
            def X_to_theta(x):
                '''
                Assumes screen distance in same units as x (a0)
                '''
                return 0.5*np.arctan2(x,screen_distance)
            # We must determine q at each point for finding I
            def X_to_q(x):
                '''
                Returns q [1/a0]
                '''
                theta = X_to_theta(x)
                lamb = E_to_lamb(self.photon_energy)
                return 4*np.pi*np.abs(np.sin(theta))/lamb
            # To get a geometry-independent plot, we may want to switch from x coords to the q component parallel to the screen.
            # We can determine this q_scr component since we know the photon momentum and theta            
            def X_to_q_scr(x):
                ''' 
                Returns component of q parallel to screen [1/a0]
                '''
                theta = X_to_theta(x)
                k0 = self.photon_momentum
                return k0*np.sin(2*theta)

            ## Calculate q for each cell.
            max_q = X_to_q(np.sqrt(2*screen_width**2))
            result.cell_width = screen_width/len(cell)
            q_grid = np.empty(cell.shape)
            result.xy = np.empty(cell.shape+(2,))
            result.I = np.zeros(cell.shape)            
            for i, x in enumerate(cell):
                for j, y in enumerate(x):
                    x = result.cell_width*(i-(len(cell)-1)/2)
                    y = result.cell_width*(j-(len(cell)-1)/2)
                    if self.hemisphere_screen:
                        print ("hemisphere not supported for SPI yet")
                    
                    dist = np.sqrt(x**2+y**2)
                    q_grid[i,j] = X_to_q(dist)
                    result.xy[i,j] = np.array([x,y])
            result.q_scr_xy = X_to_q_scr(result.xy)

            # Calculate I for each cell from q (need to vectorise q)
            for rot_x in range(self.x_orientations):
                self.x_rot_matrix = rotaxis2m(self.x_rotation,Bio_Vect(1, 0, 0))      
                #self.y_rotation = 0                 
                for rot_y in range(self.y_orientations):
                    print("Imaging at x, y, rotations:",self.x_rotation,self.y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.y_rotation,Bio_Vect(0, 1, 0))      
                    self.y_rotation += 2*np.pi/self.y_orientations              
                    result.I += self.generate_cell_intensity(q_grid,result.xy)
                    #mask = (q_grid < self.min_q) + (q_grid > self.max_q)
                    #mask = np.zeros(q_grid.shape,dtype=bool)
                    #mask[q_grid < self.min_q] = True; mask[q_grid > self.max_q] = True

                    mask = q_grid
                    mask[mask < self.min_q] = 0; mask[mask > self.max_q] = 0
                    mask[mask != 0] = 1

                    #print("MASK",mask)
                    #result.I[mask == 0] = None
                    result.I*=mask
                self.x_rotation += 2*np.pi/self.x_orientations             
            # convert to angstrom
            result.xy /= ang_per_bohr
            result.cell_width /= ang_per_bohr
            result.q_scr_xy *= ang_per_bohr
            # Average out intensity # NOTE intensities aren't aligned. Shouldn't be using R factor directly on result from multiple orientations for SPI.
            # TODO Intensity should really be a list of results.
            result.I /= (self.y_orientations*self.x_orientations)
            #result.package_up()
            return result        

        else:
            # Iterate through each orientation of crystal, picklin' up a file for each orientation
            used_orientations = []
            for j, cardan_angles in enumerate(self.orientation_set):               
                print("Imaging orientation",j)
                bragg_points, miller_indices,cardan_angles = self.bragg_points(target,cell_packing = target.cell_packing, cardan_angles = cardan_angles,random_orientation=random_orientation)
                used_orientations.append(cardan_angles)
                num_points = int(len(bragg_points))
                result = Results(num_points,j)
                # Get the q vectors where non-zero
                i = 0
                #TODO vectorise
                # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                point = self.generate_point(bragg_points,cardan_angles)
                # fix this filling stuff
                result.phi = point.phi
                result.phi_aligned = point.phi_crystal_aligned
                result.q = point.q_parr_screen
                result.I += point.I
                
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
    class Cell(Feature):
        def __init__(self,*args):
            super().__init__(*args)
            
        
    def generate_cell_intensity(self,q,xy_grid):
        solid_angle = 1  # approximate all as same for now.
        theta = self.q_to_theta(q)
        q_parr_screen = self.q_to_q_scr(q)     
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")           
        cell = self.Cell(q,q_parr_screen,theta)
        phis = np.arctan2(xy_grid[:,:,1],xy_grid[:,:,0])
        #print("phis",phis)
        return self.illuminate(cell,phis=phis,SPI_proj_solid_angle = solid_angle)

    def generate_ring(self,q,phi_array,angle_subtended):
        '''Returns the intensity(phi) array and the radius for given q.'''
        #print("q=",q)
        #r = self.q_to_r(q)# TODO
        theta = self.q_to_theta(q)
        q_parr_screen = self.q_to_q_scr(q)
        solid_angle = (2*np.pi/len(phi_array)) * (np.pi/angle_subtended)   # solid angle in fractions [sp]. Solved for sphere but will be same when project to flat detector. 
        proj_solid_angle = solid_angle  # laser source.
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")
        ring = self.Ring(q,q_parr_screen,theta)
        ring.I = self.illuminate(ring,phis=phi_array,SPI_proj_solid_angle=proj_solid_angle)
        return ring 
    def generate_point(self,G,cardan_angles): # G = vector
        ''' We store variables like G since they correspond to the special bragg points, unlike points from SPI that are just samples of the continuous pattern. 
        
        '''
        if len(G.shape) > 1:
            G = np.moveaxis(G,0,len(G.shape)-1)  # moves the axis corresponding to the individual momentum to the back, giving us dim = [3,num_G].
        q = np.sqrt(np.power(G[0],2)+np.power(G[1],2)+np.power(G[2],2))
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
        
        point.I = self.illuminate(point,cardan_angles=cardan_angles)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if q_parr_screen != check:
        #     print("Error, q_parr_screen =",q_parr_screen,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
    def illuminate(self,feature, phis = None,cardan_angles = None, SPI_proj_solid_angle=None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet.
        phis is an ndarray containing the angle of each point to calculate for (SPI) rings, but should contain only 1 element for points.
        """
        # if phis == None:
        #     phis = self.phi_array
        
        if type(feature) == self.Spot:
            SPI = False
        elif type(feature) in [self.Cell,self.Ring]:
            SPI = True

        times_used = None
        for species in self.target.species_dict.values():
            species.set_scalar_form_factor()
            times_used = species.times_used

        F_shape = tuple()
        if type(feature) is self.Spot:
            F_shape += (self.t_fineness,)           
            if len(feature.G.shape) > 1:
                F_shape += (len(feature.G[2]),) #[times,num_G]   
        else:
            if type(feature) is self.Ring:
                F_shape += phis.shape  
            F_shape += (self.t_fineness,)# [?phis?,times]
            if type(feature.q) is np.ndarray:
                F_shape += feature.q.shape          # [?phis?,times,feature.q.shape]
                
        # Technically sum of F(t)*sqrt(J(t)), where F = sum(f(q,t)*T(q)), and J(t) is the incident intensity, thus accounting for the pulse profile.
        F_sum = np.zeros(F_shape,dtype="complex_")  
        for species in self.target.species_dict.values():
            print("------------------------------------------------------------")
            print("Getting contribution to integrand from species",species.name)
            if not np.array_equal(times_used,species.times_used):
                raise Exception("Times used don't match between species.")        
            # iterate through every atom including in each symmetry of unit cell (each asymmetric unit)
            max_atoms_per_loop = 200  # Restrict array size to prevent computer tapping out.
            for s in range(len(self.target.sym_rotations)):
                print("Working through symmetry",s)
                num_atom_batches = int(len(species.coords)/max_atoms_per_loop)+1
                for a_batch in range(num_atom_batches):
                    atm_idx = np.arange(len(species.coords)*s+max_atoms_per_loop*a_batch, len(species.coords)*s + min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                    print("atoms",atm_idx[0],"-",atm_idx[-1])
                    relative_atm_idx = np.arange(max_atoms_per_loop*a_batch, min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                    # Rotate to target's current orientation 
                    # rot matrices are from bio python and are LEFT multiplying. TODO should be consistent replace this with right mult. 
                    R = np.array(species.coords[relative_atm_idx[0]:relative_atm_idx[-1]+1]) 
                    coord = self.target.get_sym_xfmed_point(R,s)  + species.error[atm_idx] # dim = [ N, 3]
                    coord = coord @ self.x_rot_matrix  
                    coord = coord @ self.y_rot_matrix
                    coord = coord @ self.z_rot_matrix
                    # Get spatial factor T
                    if SPI:
                        T = self.SPI_interference_factor(phis,coord,feature)  #[num_G]
                    else:
                        T= self.interference_factor(coord,feature,cardan_angles)   # if grid: [phis,qx,qy]  if ring: [phis,q] (unimplemented)
                    f = species.get_stochastic_f(atm_idx, feature.q)  / np.sqrt(self.target.num_cells) # Dividing by np.sqrt(self.num_cells) so that fluence is same regardless of num cells. 
                    #print(F_sum.shape,T.shape,f.shape) 
                    if SPI: 
                        if type(feature) is self.Cell:
                            F_sum += np.sum(T[:,None,...]*f,axis=0)          #[num_atoms,None,qX,qY] X [num_atoms,times,qX,qY]
                        if type(feature) is self.Ring:           
                            F_sum += np.sum(T[:,None] * f[None,...],axis=0)        # [?phis?,momenta,None]X[None,times, feature.q.shape]  -> [?phis?,times,feature.q.shape] if q is an array
                    else:
                        F_sum += np.sum(T[None,:] * f)                           # [None,num_G]X[times,num_G]  ->[times,num_G] 
                    #print("s,F_sum",s,F_sum)
        if (times_used[-1] == times_used[0]):   
            raise Exception("Intensity array's final time equals its initial time") #I =  np.square(np.abs(F_sum[:,0]))  # 
        else:
            time_axis = 0
            if type(feature) is self.Ring:
                time_axis = 1
            I = np.trapz(np.square(np.abs(F_sum)),times_used,axis = time_axis) / (times_used[-1]-times_used[0])       #[num_G] for points, or for SPI: [phis,feature.q.shape], corresponding to rings or square grid
        
        # For unpolarised light (as used by Neutze 2000). (1/2)r_e^2(1+cos^2(theta)) is thomson scattering - recovering the correct equation for a lone electron, where |f|^2 = 1 by definition.    
        thet = self.q_to_theta(feature.q)
        r_e_sqr = 2.83570628e-9        
        I*= r_e_sqr*(1/2)*(1+np.square(np.cos(2*thet)))
        
        print("Total screen-incident intensity = ","{:e}".format(np.sum(I)))

        if SPI:
            # For crystal we can approximate infinitely small pixels and just consider bragg points. (i.e. The intensity will be the same so long as the pixel isn't covering multiple points.)
            # But for SPI need to take the pixel's size into account. Neutze 2000 makes the following approximation:
            I *=  SPI_proj_solid_angle    # equiv. to *= solid_angle
        #print(I)
        return I

    def illuminate_average(self,feature,phis = None,cardan_angles = None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet."""
        if phis == None:
            phis = self.phi_array
        F = np.zeros(phis.shape,dtype="complex_")
        for species in self.target.species_dict.values():
            species.set_scalar_form_factor(feature.q)
            # iterate through each symmetry of unit cell (each asymmetric unit)
            for s in range(len(self.target.sym_rotations)):
                for R in species.coords:
                        # Rotate to crystal's current orientation 
                        R = R.left_multiply(self.y_rot_matrix)  
                        R = R.left_multiply(self.x_rot_matrix)   
                        # PDB format note: if x axis has dim of X, we have a point between [- X, X].
                        coord = np.multiply(R.get_array(),self.target.sym_rotations[s]) + np.multiply(self.target.cell_dim,self.target.sym_translations[s])
                        # Get spatial factor T
                        T = np.zeros(phis.shape,dtype="complex_")
                        if SPI:
                            T = self.SPI_interference_factor(phis,coord,feature)
                        else:
                            T= self.interference_factor(coord,feature,cardan_angles)
                        F += species.ff_average*T   
                        # Rotate atom for next sample            
        I = np.square(np.abs(F))        
        

        return I 

    def interference_factor(self,coord,feature,cardan_angles,axis=0):
        """ theta = scattering angle relative to z-y plane """ 
        # Rotate our G vector BACK to the real laser orientation relative to the crystal.
        q_vect = self.rotate_G_to_orientation(feature.G.copy(),*cardan_angles,inverse=True)[0]       
        q_dot_r = np.apply_along_axis(np.dot,axis,q_vect,coord) # dim = [num_G]                      
        T = np.exp(-1j*q_dot_r)
        return T

    def SPI_interference_factor(self,phi_array,coord,feature):  #TODO refactor SPI features so can just use above (allows for realignment)
        """ theta = scattering angle relative to z-y plane 
        
        """ 
        q_z = np.multiply(feature.q,np.sin(feature.theta))                          # dim = [qX,qY] for square screen with cell coords given by X,Y 
        # phi is angle of vector relative to x-y plane, pointing from screen centre to point of incidence. 
        # TODO need to get working with rings again since they were designed to have multi phi for each q. Here it is 1:1
        if type(feature) is not self.Cell:
            raise Exception("accidentally removed ring implementation sorry!")
        q_z = q_z       
        q_y = feature.q_parr_screen * np.sin(phi_array)
        q_x = feature.q_parr_screen * np.cos(phi_array)   # dim = [qX,qY]
        #print("====="); print(phi_array);print("-----");print(q_z.shape,q_y.shape,q_x.shape)   
        q_vect = np.array([q_x,q_y,q_z])
        q_vect = np.moveaxis(q_vect,0,len(q_vect.shape)-1)    # dim = [qX,qX,3]
        coord = np.moveaxis(coord,0,-1)
        q_dot_r = np.apply_along_axis(np.matmul,len(q_vect.shape)-1,q_vect,coord) # dim = [phis,qX,qY]       q_dot_r = np.tensordot(q_vect,coord,axes=len(q_vect.shape)-1)
        q_dot_r = np.moveaxis(q_dot_r,-1,0)
        coord = np.moveaxis(coord,-1,0)
        if (type(feature.q) == np.ndarray): 
            if len(coord.shape) == 2 and q_dot_r.shape != (coord.shape[0],) + feature.q.shape or len(coord.shape) == 1 and q_dot_r.shape != feature.q.shape:
                raise Exception("Unexpected q_dot_r shape.",q_dot_r.shape)
        T = np.exp(-1j*q_dot_r) 
        return T   
    def bragg_points(self,crystal, cell_packing, cardan_angles,random_orientation=False):
        ''' 
        Using the unit cell structure, find non-zero values of q for which bragg 
        points appear.
        lattice_vectors e.g. = [a,b,c] - length of each spatial vector in orthogonal basis.
        Currently only checked to work with cubics. 
        '''

        def get_G(miller_indices,cartesian=True):
            '''
            Get G in global cartesian coordinates (crystal/zone axes not implemented),
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
                G, used_angles = self.rotate_G_to_orientation(G,*cardan_angles,random = random_orientation,G_idx_first = True) 
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
        q_1_max = self.max_q
        q_2_max = self.max_q
        q_3_max = self.max_q        # q = (0,0,l) case.
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
        # If we used a random orientation, we now lock in the orientations just generated and contained in cardan_angles.
        random_orientation = False 

        if self.max_triple_miller_idx != None:
            m = self.max_triple_miller_idx
            max_g_vect = get_G(np.full((1,3),m))[0][0]
            self.max_q = min(self.max_q,np.sqrt(((max_g_vect[0])**2+(max_g_vect[1])**2+(max_g_vect[2])**2)))
        
        print("using q range of ", self.min_q*1.88973,"-",self.max_q*1.88973," angstrom-1")
        min_max_q_rule = lambda g: self.min_q <= np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2)) <= self.max_q
        #max_q_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.max_q
        q0 = self.photon_momentum
        
        # Catch the miller indices with a boolean mask
        mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(min_max_q_rule,1,G_temp)*np.apply_along_axis(self.mosaic_elastic_condition,1,G_temp)
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

    def rotate_G_to_orientation(self,G,alpha=0,beta=0,gamma=0,random=False,inverse = False, G_idx_first = False):
        '''
        Assumes G.shape = (num_momenta,3)
        Rotates the G vectors to correspond to the crystal when oriented to the given cardan angles. 
        alpha: x-axis angle [radians].
        beta:  y-axis angle [radians].
        gamma: z-axis angle [radians].
        ## Returns:
        G: transformed G
        list: [alpha,beta,gamma] (the angles used)
        '''
        moved_axis = False
        if G_idx_first and G.shape[len(G.shape)-1] != 3:
            raise Exception("please give array of form [...,num_G,3]")
        elif not G_idx_first and G.shape[len(G.shape) -2] !=3:
            raise Exception("please give array of form [...,3,num_G]")
            
        if G_idx_first:
            G = np.swapaxes(G,len(G.shape)-2,len(G.shape)-1)

        if type(random) != bool:
            raise Exception("random_orientation must be boolean")
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
        G = np.around(rot_matrix,decimals=10) @ G     # [3,3]X[3,num_G]

        if G_idx_first:
            G = np.swapaxes(G,len(G.shape)-2,len(G.shape)-1)  
        return G, [alpha,beta,gamma]  # [num_G,3]
                
    def q_to_theta(self,q):
        '''
        q momentum transfer [1/bohr]
        '''
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

    #TODO I'll have to check later if this is crazy or not
    '''   
    def q_to_scr_pos(self,q):
        """ 
        """
        theta = self.q_to_theta(q)
        # From here it will help to consider q as referring to the final momentum. /// what is this comment??????? 
        v_x = c_au*np.cos(2*theta)
        v_y = c_au*np.sin(2*theta)
        D = self.detector_distance
        t = D/v_x
        radius = v_y*t
        return radius
    '''
    def q_to_q_scr(self,q):
            """ Returns the screen-parallel component of q (or equivalently the final momentum).
            """
            theta = self.q_to_theta(q)
            return q/(np.sin(np.pi/2-theta))   
    def q_to_q_scr_curved(self,G):
        return np.sqrt(G[0]**2+G[1]**2)
    
    def mosaic_elastic_condition(self,q_vect):
        ''' determines whether vector q is allowed.'''
        # Rocking angle version
        q = np.sqrt(q_vect[0]**2 + q_vect[1]**2 + q_vect[2]**2)
        if q > self.max_q:
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

def scatter_scatter_plot(neutze_R = True, crystal_aligned_frame = False ,SPI_result1 = None, SPI_result2 = None, full_range = True,num_arcs = 50,num_subdivisions = 40, result_handle = None, compare_handle = None, normalise_intensity_map = False, show_grid = False, cmap_power = 1, cmap = None, min_alpha = 0.05, max_alpha = 1, bg_colour = "black",solid_colour = "white", show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False,  fixed_dot_size = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None):
    ''' (Complete spaghetti at this point.)
    Plots the simulated scattering image.
    result_handle:

    compare_handle:
        results/compare_handle/ is the directory of results that will be subtracted from those in the results directory.
        Must have same orientations.
    
    '''
    print("=====================Plotting===========================")
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
        #plt.gcf().set_figwidth(20)         # if not using widget magic
        #plt.gcf().set_figheight(20)        #       
            
        plt.gcf().set_figwidth(10)
        plt.gcf().set_figheight(10)                  
    
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

        if compare_handle != None:
            if log_I:
                print("Not plotting logarithmic I, not supported for comparisons.")
                log_I = False  
            
            all_I_ideal = []
            all_I_real = []
        
        # Iterate through each orientation (file) to get minimum/maximum for normalisation of plot of scattering image (not difference image).
        max_z = -np.inf
        min_z = np.inf

        fig = plt.figure()
        # fig.canvas.layout.width = '40%'
        # fig.canvas.layout.height = '40%'
        # fig.canvas       
        ax = fig.add_subplot(projection="polar")        
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
                radial_axis = result1.q        
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
                if (radial_axis[i], phi[i],result1.image_index)  in processed_copies:
                    unique_values_mask[i] = False
                    continue
                else:
                    unique_values_mask[i] = True
                # if (radial_axis[i] > radial_lim):
                #     unique_values_mask[i] = False
                    
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
                    all_I_real.extend(I1)
                    all_I_ideal.extend(I2)


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
                print("Neutze R:")
                sqrt_ideal = np.sqrt(all_I_ideal)
                sqrt_real = np.sqrt(all_I_real)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)
                print("Regular R:")
                R = np.sum(np.abs((sqrt_real - sqrt_ideal)))/np.sum(sqrt_ideal)
                print(R)
                print("sum of real","{:e}".format(np.sum(sqrt_real)),"sum of ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(np.sum(np.abs((sqrt_real - sqrt_ideal)))))
                #neutze_histogram = np.histogram2d(phi, radial_axis, weights=R, bins=(np.array([-np.pi,np.pi]), radial_edges))[0]             
                #plot_sectors(neutze_histogram)
    
    ## Continuous (SPI)
    else:
        result1 = SPI_result1
        result2 = SPI_result2        
        ## Square grid
        if len(result1.I.shape) > 1:#if type(result1) == Results_Grid:
            if log_I: 
                z1 = np.log(result1.I)
                z1[result1.I == 0] = None
                if result2 != None:
                    z2 = np.log(result2.I)
                    z2[result2.I == 0] = None
            else:
                #TODO fix up this masked stuff i didnt finish 
                # z1 -= cutoff_log_intensity
                # z1[z1<0] = 0
                z1 = result1.I.copy()        
                if result2 != None:
                    #z2 -= cutoff_log_intensity
                    #z2[z2<0] = 0
                    z2 = result2.I.copy()        
            if log_I and cutoff_log_intensity != None:
                np.ma.array(z1, mask=(z1<cutoff_log_intensity)*(np.isnan(z2)))
                if result2 != None:
                    z2 = np.ma.array(z2, mask=(z2<cutoff_log_intensity)*(np.isnan(z2)))
            else:
                z1 = np.ma.array(z1,mask = np.isnan(z1))
                z2 = np.ma.array(z2,mask = np.isnan(z2)) 
            print(z1)
            print(z2)       
            print("Result 1 (Damaged):")
            print("Total screen-incident intensity:","{:e}".format(np.sum(result1.I)))
            if result2 != None:
                combined_data = np.array([z1,z2])
            else: 
                combined_data = z1
            z_min, z_max = np.nanmin(combined_data), np.nanmax(combined_data)       
            current_cmap = plt.colormaps.get_cmap("viridis")
            current_cmap.set_bad(color='black')
            z1_map = plt.imshow(z1,vmin=z_min,vmax=z_max,cmap=current_cmap)
            plt.colorbar(z1_map)
            plt.show()
            if result2 != None:
                print("Result 2 (Undamaged):")
                print("Total screen-incident intensity:","{:e}".format(np.sum(result2.I)))
                z2_map = plt.imshow(z2,vmin=z_min,vmax=z_max,cmap=current_cmap)
                plt.colorbar(z2_map)
                plt.show()
                print("Neutze R: ")
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R_cells = np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal))
                R = np.sum(R_cells)
                print(R) 
                R_cells *= len(R_cells)**2# multiply by num cells (such that R is now the average) 
                R_map = plt.imshow(R_cells,cmap=cmap)     
                plt.colorbar(R_map)
                ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_scr_xy[0]]
                plt.xticks(ticks,ticklabels)
                plt.yticks(ticks,ticklabels)
                plt.show()

                print("Regular R:")
                A = np.abs(sqrt_real - sqrt_ideal)
                R_num = np.sum(A)
                R_den = np.sum(sqrt_ideal)
                R = R_num/R_den                
                print(R)   
                print("sum of real","{:e}".format(np.sum(sqrt_real)),"sum of ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(R_num))       
        ## Circle grid
        else:
            if log_I: 
                z = np.log(result1.I)
            else:
                z = result1.I        
            if log_I and cutoff_log_intensity != None:
                #cutoff_log_intensity = -1
                z -= cutoff_log_intensity
                z[z<0] = 0
                pass
            #radial_axis = result.r
            if plot_against_q:
                radial_axis = result1.q         
            fig = plt.figure()
            ax = fig.add_subplot(projection="polar")
            # phi_mesh = result.phi_mesh
            # if crystal_aligned_frame:
            #     phi_mesh = result.phi_aligned_mesh   
            #     print("crystal aligned not implemented for SPI yet")    
            print(result1.phi.shape)
            print(radial_axis.shape)
            print(z.shape) 
            ax.pcolormesh(result1.phi, radial_axis, z,cmap=cmap)
            ax.plot(result1.phi, radial_axis, color = 'k',ls='none')
            plt.grid(alpha=min(show_grid,0.6),dashes=(5,10)) 

            if radial_lim:
                bottom,top = plt.ylim()
                plt.ylim(bottom,radial_lim)
            if log_radial:
                plt.yscale("log")

            if result2 != None:
                print("Neutze R:")
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)            
                print("Regular R:")
                R = np.sum(np.abs((sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)

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




##### https://scripts.iucr.org/cgi-bin/paper?S0021889807029238, http://superflip.fzu.cz/

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

#####
# stylin' 
def stylin(SPI=False,SPI_max_q=None,SPI_result1=None,SPI_result2=None):
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
    try:
        cmap = shiftedColorMap(matplotlib.cm.RdYlGn_r,midpoint=0.2)#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
    except: 
        cmap =  plt.get_cmap("shiftedcmap")
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
    q_scr_max = experiment1.max_q#experiment1.max_q*(1/np.sqrt(2)) # for flat screen. as q_scr = qcos(theta). (q_z = qsin(theta) =ksin^2(theta)), max theta is 45. (Though experimentally ~ 22 as of HR paper)

    zoom_to_fit = False
    if not zoom_to_fit:
        if use_q:
            radial_lim = q_scr_max+0.2
        #els:
    #else:
    #     radial_lim = screen_radius#min(screen_radius,experiment.q_to_r(experiment.max_q))
    #     print(radial_lim)
    #     if use_q:
    #         #radial_lim = experiment.r_to_q_scr(radial_lim)
    #         radial_lim = q_scr_lim #radial_lim = min(q_scr_lim,experiment.q_to_q_scr(experiment.max_q))
    #         print(radial_lim)
    # else:
    #     radial_lim = None

    #TODO fix above to work with distance

    # R Sectors
    if not SPI:
        scatter_scatter_plot(crystal_aligned_frame = False,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
        # Intensity of experiment 1. 
        scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
        # R Sectors aligned
        scatter_scatter_plot(crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
        # Intensity of experiment 1 aligned. 
        scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)
    else:
        #TODO Get the colors to match neutze.
        use_q = True
        log_radial = False
        log_I = True
        cutoff_log_intensity = None # -1 or None
        cmap = "PiYG_r"# "plasma"
        radial_lim = SPI_max_q
        scatter_scatter_plot(log_I = log_I, cutoff_log_intensity = cutoff_log_intensity, SPI_result1=SPI_result1,SPI_result2=SPI_result2,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap)

    #NEED TO CHECK. We have a 1:1 mapping from q to q_parr, but with our miller indices we are generating multiple q_parr with diff q.
    # So we SHOULD get the same q_parr with different q_z. Which makes sense since we are just doing cosine. But still icky maybe?
    # Need to double check we get different intensities for same q_parr. Pretty sure that's implemented.
#stylin()

def res_to_q(d):
    '''
    Pass d = None to use default resolution (when put into exp_args)
    '''
    if d == None:
        return None
    return 2*np.pi/d

#TODO 
# log doesnt work atm.
# Get the rings to correspond to actual rings

##%% vvvvvvvvv
### Simulate
#TODO
# have max q, that matches neutze
# implement rhombic miller indices as the angle is actually 120 degrees on one unit cell lattice vector (or just do SPI)


#  #TODO make this nicer and pop in a function 
DEBUG = False

target_options = ["neutze","hen","tetra"]
#============------------User params---------==========#

target = "tetra" #target_options[2]
top_resolution = 2
bottom_resolution = None#30

#### Individual experiment arguments 
start_time = -6
end_time = 2.8#10 #0#-9.80    
laser_firing_qwargs = dict(
    SPI = True,
    pixels_across = 100,  # shld go on xfel params.
)
##### Crystal params
crystal_qwargs = dict(
    cell_scale = 2,  # for SC: cell_scale^3 unit cells 
    positional_stdv = 0,#0.2, # RMS in atomic coord position [angstrom]
    include_symmetries = True,  # should unit cell contain symmetries?
    cell_packing = "SC",
    rocking_angle = 0.3,  # (approximating mosaicity)
    orbitals_as_shells = True,
    #CNO_to_N = True,   # whether the laser simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
)

#### XFEL params
tag = "" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
#TODO make it so reflections don't overwrite same orientation, as stochastic now.
random_orientation = True # if true, overrides orientation_set with random orientations (only affects non-SPI)
energy = 6000 # eV
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

#orientations
num_orients = 5
# [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Overridden if random orientations.
ax_x = 1
ax_y = 1
ax_z = 0


# Optional: Choose previous folder for crystal results
chosen_root_handle = None # None for new. Else use "tetra_v1", if want to add images under same params to same results.
#=========================-------------------------===========================#

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
    
#----------orientation thing that needs to be to XFEL class (TODO)-----------------------#

if ax_x == ax_y == ax_z and ax_z == 0 and random_orientation == False:
    throwabxzd
# if not random_orientation:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print("ori set if not random:", orientation_set)  #TODO double check shows when random and if so disable when random
#---------------------------------#
if target == "neutze": #T4 virus lys
    pdb_path = "/home/speno/AC4DC/scripts/scattering/2lzm.pdb"
    cell_dim = np.array((61.200 ,  61.200 ,  61.2))#TODO should read this automatically...
    target_handle = "D_lys_neutze_simple_7"  
    allowed_atoms_1 = ["N_fast","S_fast"]
    #allowed_atoms_2 = ["N_fast","S_fast"]      
    CNO_to_N = True
elif target == "hen": # egg white lys
    pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb"
    cell_dim = np.array((79.000  , 79.000  , 38.000))
    target_handle = "D_lys_neutze_simple_7"  
    allowed_atoms_1 = ["N_fast","S_fast"]
    #allowed_atoms_2 = ["N_fast","S_fast"]      
    CNO_to_N = True
elif target == "tetra": 
    pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb" 
    cell_dim = np.array((4.813 ,  17.151 ,  29.564))
    target_handle = "carbon_gauss_32"
    allowed_atoms_1 = ["C_fast"]
    #allowed_atoms_2 = ["C_fast"]      
    CNO_to_N = False
else:
    raise Exception("'target' invalid")
#-------------------------------#

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
    SPI_result1 = experiment1.spooky_laser(start_time,end_time,target_handle,crystal, random_orientation = random_orientation, **laser_firing_qwargs)
    SPI_result2 = experiment2.spooky_laser(start_time,end_time,target_handle,crystal_undmged, random_orientation = False, **laser_firing_qwargs)
    stylin(SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
else:
    exp1_orientations = experiment1.spooky_laser(start_time,end_time,target_handle,crystal, random_orientation = random_orientation, **laser_firing_qwargs)
    create_reflection_file(exp_name1)
    experiment2.orientation_set = exp1_orientations  # pass in orientations to next sim, random_orientation must be false so not overridden!
    experiment2.spooky_laser(start_time,end_time,target_handle,crystal_undmged, random_orientation = False, **laser_firing_qwargs)
    stylin()

#%%
orb_occs1, times_used = crystal.ff_calculator.random_state_snapshots("C_fast") 
orb_occs2, times_used = crystal.ff_calculator.random_state_snapshots("C_fast") 
twoatoms = np.array(orb_occs1,orb_occs2)
crystal.ff_calculator.random_states_to_f_snapshots(times_used,twoatoms, 1, "C_fast")     
#%%
if laser_firing_qwargs["SPI"]:
    stylin(SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
else:
    stylin()
#^^^^^^^
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
### Simulate crambin
tag += 1
#  #TODO make this nicer and pop in a function 
DEBUG = False
energy =17445   #crambin - from resolutions. in pdb file, need to double check calcs: #q_min=0.11,q_cmaxres = 3.9,  my defaults: pixels_per_ring = 400, num_rings = 200


exp_name1 = str(root) + "1_" + str(tag)
exp_name2 = str(root) + "2_" + str(tag)



allowed_atoms_1 = ["C_fast","N_fast","O_fast","S_fast"]
allowed_atoms_2 = ["C_fast","N_fast","O_fast","S_fast"]

start_time = -10
end_time = 10#-9.95

#orientation_set = [[5.532278012665244, 1.6378991611963682, 1.3552062099726534]] # Spooky spider orientation
num_orients = 30
# [ax_x,ax_y,ax_z] = vector parallel to axis. Overridden if random orientations.
ax_x = 0
ax_y = 0
ax_z = 1
random_orientation = True  # if true, overrides orientation_set

rock_angle = 0.15 # degrees
screen_type = "hemisphere"

pdb_path = "/home/speno/AC4DC/scripts/scattering/3u7t.pdb" #Crambin

target_handle = "Improved_Lys_mid_6"

q_cutoff = None

test_the_onion = False  # <- a test case I'm familiar with, I'm using it to check if anything breaks.
#---------------------------------#


# if not random:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

if test_the_onion:
    orientation_set = [[0,0,0],[np.pi/4,0,0]]
    #max_q = 10

experiment1 = XFEL(exp_name1,energy,100, q_cutoff = q_cutoff,screen_type = screen_type, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)
experiment2 = XFEL(exp_name2,energy,100, q_cutoff = q_cutoff,screen_type = screen_type, orientation_set = orientation_set, pixels_per_ring = 400, num_rings = 400,t_fineness=100)


crystal_dmged = Crystal(pdb_path,allowed_atoms_1,rocking_angle = rock_angle*np.pi/180,CNO_to_N=False,cell_packing = "SC")
crystal_dmged.set_cell_dim(22.795 ,  18.826 ,  41.042)
if test_the_onion:
   crystal_dmged.set_cell_dim(20,20,20)
crystal_dmged.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0]))   #2555

crystal_undmged = copy.deepcopy(crystal_dmged)
crystal_undmged.is_damaged = False


SPI = False

exp1_orientations = experiment1.spooky_laser(-10,end_time_1,target_handle,crystal, random_orientation = random_orientation, SPI=SPI)
create_reflection_file(exp_name1)
experiment2.orientation_set = exp1_orientations
experiment2.spooky_laser(-10,end_time_2,target_handle,crystal, random_orientation = False, SPI=SPI)

stylin()




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

target_handle = "Improved_Lys_mid_6"

#---------------------------------#


# if not random:
orientation_axis = Bio_Vect(ax_x,ax_y,ax_z)
orientation_set = [rotaxis2m(angle, orientation_axis) for angle in np.linspace(0,2*np.pi,num_orients,endpoint=False)]
orientation_set = [Rotation.from_matrix(m).as_euler("xyz") for m in orientation_set]
print(orientation_set)

experiment = XFEL("SPI",energy,100, screen_type = "hemisphere_screen", pixels_per_ring = 500, num_rings = 500,t_fineness=100)
crystal = Crystal(pdb_path,allowed_atoms,CNO_to_N=False,cell_packing = "SC")

SPI = True

SPI_result = experiment.spooky_laser(-10,end_time,target_handle,crystal, SPI=SPI)

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
stylin()
# %%
