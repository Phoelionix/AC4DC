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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from plotter_core import Plotter

class Results():
    pass

class XFEL():
    def __init__(self, photon_energy, detector_distance,  x_orientations = 1, y_orientations = 1, q_min=0.5,q_max=2.4, pixels_per_ring = 400, num_rings = 50,t_fineness=100):
        """ #### Initialise the imaging experiment's controlled parameters
        photon_energy [eV]:
            Should be the same as that given in the original input file!
        detector_distance [bohr];
            The distance between the target and the centre of the detector 
        q_max [1/bohr]:
        pixels_per_ring: 
            The number of different values for alpha to plot. Alpha is the angle that the incident light makes with 
            the y axis, (z-axis is firing line. y-axis is any arbitrary axis perpendicular to z-axis.)
        num_rings: 
            Determines the number of q to calculate. Note pixels are currently just points
        t_fineness:
            Number of time steps to calculate.
        """
        self.photon_energy = photon_energy
        self.detector_distance = detector_distance
        self.q_min = q_min
        self.q_max = q_max
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.y_orientations = y_orientations
        self.x_orientations = x_orientations

        self.t_fineness = t_fineness

        self.alpha_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)
        self.y_rotation = 0 # Current rotation of crystal (y axis currently)
        self.x_rotation = 0 # Current rotation of crystal (y axis currently)
        self.y_rot_matrix = rotaxis2m(self.y_rotation,bio_vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.x_rotation,bio_vect(1, 0, 0))     

    def set_atomic_species(self, pl, pdb_fpath, allowed_atoms, CNO_to_N = False, orbitals_as_shells = True):
        parser=PDBParser(PERMISSIVE=1)
        structure_id = os.path.basename(pdb_fpath)
        structure = parser.get_structure(structure_id, pdb_fpath)
        model = structure[0]
        if len(structure) > 1:
            print("WARNING: Only captured first of",structure.size(),"structures!")

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
        for chain in model.get_list():
            atoms_gen = chain.get_atoms()
            for atom in atoms_gen:
                # Get ze data
                R = atom.get_vector()
                name = AC4DC_name(atom.element)
                # Pop into our desired format
                if name == None:
                    continue
                if name not in species_dict.keys():
                    species_dict[name] = self.Atomic_Species(name,pl) 
                species_dict[name].add_atom(R)
        for str in allowed_atoms:
            if str not in species_dict.keys():
                print("Warning: no",str,"atoms found.")
        self.species_dict = species_dict 


            
    class Atomic_Species():
        def __init__(self,name,pl):
            self.name = name 
            self.pl = pl 
            self.ff = 0 # form_factor
            self.coords = []  # coord of each atom in species
        def add_atom(self,vector):
            self.coords.append(vector)
        def set_scalar_form_factor(self,q):
            # F_i(q) = f(q)*T_i(q), where f = self.ff is the time-integrated average 
            self.ff = self.pl.f_average(q,self.name)  
                      
    def get_pl(self,start_time,end_time,output_handle):
        pl = Plotter(output_handle,"y")
        plt.close()
        pl.initialise_coherence_params(start_time,end_time,self.q_max,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return pl
    
    def firin_mah_lazer(self, start_time, end_time, output_handle, pdb_fpath, allowed_atoms="All", CNO_to_N = False,log=True):
        """ 
        end_time: The end time of the photon capture in femtoseconds. Not a real thing experimentally, but useful for choosing 
        a level of damage. Explicitly, it is used to determine the upper time limit for the integration of the form factor.
        pdb_fpath: The pdb file's path. Changes the variable self.atoms.
        y_orientations: 
            Number of unique y axis rotations to sample crystal. x_axis_rotations not implemented (yet?). 
        """
        pl = self.get_pl(start_time,end_time,output_handle) #TODO rename pl, maybe move method to different class.      
        self.set_atomic_species(pl,pdb_fpath,allowed_atoms,CNO_to_N)

        ring = np.empty(self.num_rings,dtype="object")

        q_samples = np.linspace(self.q_min,self.q_max,self.num_rings)

        result = Results()
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


                # Initialise stuff that is constant between images 
                if rot_y  == 0 and rot_x == 0:
                    azm = self.alpha_array
                    result.z = 0     
                    radii = np.zeros(self.num_rings)
                    for i in range(len(ring)):
                        radii[i] = ring[i].R           
                    r, alph = np.meshgrid(radii, azm)     
                    q_for_plot, alph = np.meshgrid(q_samples,azm)

                z = np.zeros(r.shape)  # z is the intensity of the plot colour.
                for ang in range(len(z)):
                    for pos in range(len(z[ang])):
                        z[ang][pos] = ring[pos].I[ang]
                        if log:
                            z[ang][pos] = np.log(z[ang][pos])                    
            
                result.z += z/(self.y_orientations*self.x_orientations)
            self.x_rotation += 2*np.pi/self.x_orientations
        result.r = r
        result.alph = alph
        result.azm = azm
        result.q = q_for_plot
        self.x_rotation = 0
        return result
    
    def plot_pattern(self,result,plot_against_q=False,log=False):
        # https://stackoverflow.com/questions/36513312/polar-heatmaps-in-python
        fig = plt.figure()
        ax = Axes3D(fig)        
        radial_axis = result.r
        if plot_against_q:
            radial_axis =  result.q
        plt.subplot(projection="polar")
        plt.pcolormesh(result.alph, radial_axis, result.z)
        plt.plot(result.azm, radial_axis , color='k', ls='none')
        if log:
            plt.yscale("log") 
        plt.grid()  # Make the grid lines represent one unit cell (when implemented).         
    
    class Ring:
        def __init__(self,q,theta):
            self.q = q
            self.theta = theta
    def generate_ring(self,q):
        '''Returns the intensity(alpha) array and the radius for given q.'''
        #print("q=",q)
        ring = self.Ring(q,self.q_to_theta(q))
        ring.R = self.q_to_x(q)
        ring.I = self.illuminate(ring)
        return ring 

    # Not crystalline yet
    def illuminate(self,ring):
        """Returns the intensity at q. Not crystalline yet."""
        F = np.zeros(self.alpha_array.shape,dtype="complex_")
        for species in self.species_dict.values():
            species.set_scalar_form_factor(ring.q)
            count = False
            for R in species.coords:
                # Rotate to crystal's current orientation 
                R = R.left_multiply(self.y_rot_matrix)  
                R = R.left_multiply(self.x_rot_matrix)   
                # Get spatial factor T
                T = np.zeros(self.alpha_array.shape,dtype="complex_")
                T= self.spatial_factor(self.alpha_array,R,ring)
                F += species.ff*T
                # Rotate atom for next sample            
        I = np.square(np.abs(F))
        

        return I            

    def spatial_factor(self,alpha_array,R,ring):
        """ theta = scattering angle relative to z-y plane """ 
        q_z = ring.q*np.sin(ring.theta) # not cos because q is hypotenuse - not perpendicular to x-y plane.
        q_z = np.tile(q_z,len(alpha_array))
        # alpha angle of vector relative to x-y plane, pointing from screen centre to point hit. 
        q_y = q_z*np.sin(alpha_array)
        q_x = q_z*np.cos(alpha_array)
        q_vect = np.column_stack([q_x,q_y,q_z])
        return np.exp(-1j*np.dot(q_vect,R.get_array()))   
                
    # E = photon energy in eV
    def q_to_theta(self,q_abs):
        E = self.photon_energy
        bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h_ev = 6.582119569e-16 *2 *np.pi  #TODO check why angstrom_per_nm is used instead by HR but with bohr units.
        lamb =  h_ev*c_nm/E
        lamb *= bohr_per_nm   
        return np.arcsin(lamb*q_abs/(4*np.pi))    

    def q_to_x(self,q_abs):
        E = self.photon_energy
        D = self.detector_distance
        theta = self.q_to_theta(q_abs)
        return D*np.tan(theta)    


    # Crystalline.
    def illuminate_miller(self, atoms):
        pass
        #non-zero q: miller indices q_hkl

        # Determine bragg peaks location for q = 4*pi*sin(theta)/lambda. 
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. (ideal) resolution.
        
        # q_max = 10
        # peaks = [0]
        # n = 0
        # while peaks[-1] < q_max:
        #     n += 1
        #     peaks.append(4*)


        # for q in peaks:
        #     q = 4*np.pi 
            

        # return intensity
    # def x_to_q(x,E,D):
    #     # D*np.tan(np.arcsin(lamb*q_abs/(4*np.pi))) = x, q = sin(arctan((x/D)))/lamb*4*np.pi
    #     bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h_ev = 6.582119569e-16 *2 *np.pi  #TODO check why angstrom_per_nm is used instead by HR but with bohr units.
    #     lamb =  h_ev*c_nm/E
    #     lamb *= bohr_per_nm       
    #     q = np.sin(np.arctan((x/D)))/lamb*4*np.pi
    #     return q

    '''
    def q_to_theta(q_abs,E):
        bohr_per_nm = 18.8973; c_nm = 2ww.99792458e17; h_ev = 6.582119569e-16 *2 *np.pi
        lamb =  h_ev*c_nm/E
        print(lamb)
        lamb *= bohr_per_nm    
        print(lamb)
        return np.arcsin(lamb*q_abs/(4*np.pi))*180/np.pi
    '''

queue = 0.4
test = XFEL(12000,100)
pl_t = test.get_pl(-10,-9.95,"Naive_Lys_C_7")
test.set_atomic_species(pl_t,"/home/speno/AC4DC/scripts/scattering/4et8.pdb",["N_fast","S_fast"],CNO_to_N=True)
print(test.generate_ring(queue).I[0])

#%% Lysozyme
experiment = XFEL(6000,100,x_orientations=1, y_orientations=1,q_max=2.4, pixels_per_ring = 400, num_rings = 750,t_fineness=100)
#%%
# Nitrogen + Sulfur
allowed_atoms_1 = ["N_fast","S_fast"]
end_time_1 = -9.95
output_handle = "Naive_Lys_C_7"
pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb"
result1 = experiment.firin_mah_lazer(-10,end_time_1,output_handle,pdb_path,allowed_atoms_1,CNO_to_N=True)
experiment.plot_pattern(result1)
#%%
allowed_atoms_2 = ["N_fast","S_fast"]
end_time_2 = -9.76
output_handle = "Naive_Lys_C_7"
pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb"
result2 = experiment.firin_mah_lazer(-10,end_time_2,output_handle,pdb_path,allowed_atoms_2,CNO_to_N=True)
experiment.plot_pattern(result2)

#%% Difference
result3 = Results()
result3.r = result1.r
result3.q = result1.q
result3.z = result1.z-result2.z
result3.alph = result1.alph
result3.azm = result1.azm
experiment.plot_pattern(result3)
#
#%% Tetrapeptide 
experiment = XFEL(6000,100,x_orientations = 1, y_orientations=1,q_min=0.2,q_max=2.4, pixels_per_ring = 200, num_rings = 250,t_fineness=100)
use_q = True
use_log = False
pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb"

# 1
allowed_atoms_1 = ["C_fast"]
end_time_1 = -9.99
output_handle = "B_tetrapeptide_1"
result1 = experiment.firin_mah_lazer(-10,end_time_1,output_handle,pdb_path,allowed_atoms_1,CNO_to_N=False)
experiment.plot_pattern(result1,plot_against_q = use_q,log=use_log)
#%% 2
allowed_atoms_2 = ["C_fast"]
end_time_2 = -5
output_handle = "C_tetrapeptide_2"
result2 = experiment.firin_mah_lazer(-10,end_time_2,output_handle,pdb_path,allowed_atoms_2,CNO_to_N=False)
experiment.plot_pattern(result2,plot_against_q = use_q,log=use_log)
#%% Difference
result3 = Results()
result3.r = result1.r
result3.q = result1.q
result3.z = np.abs(np.sqrt(np.exp(result1.z))-np.sqrt(np.exp(result2.z)))/np.sqrt(np.exp(result1.z))
result3.alph = result1.alph
result3.azm = result1.azm
experiment.plot_pattern(result3,plot_against_q = use_q,log=use_log)


#TODO 
# log doesnt work when atm.
# %%
