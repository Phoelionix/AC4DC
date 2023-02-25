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
from plotter_core import Plotter

c_au = 137.036; eV_per_Ha = 27.211385 

class Results():
    pass


class Crystal():
    def __init__(self, pdb_fpath, allowed_atoms, CNO_to_N = False, orbitals_as_shells = True):
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
        self.cell_dim = np.array([x,y,z])
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
                      
    def get_ff_calculator(self,start_time,end_time,damage_output_handle):
        ff_calculator = Plotter(damage_output_handle,"y")
        plt.close()
        ff_calculator.initialise_coherence_params(start_time,end_time,self.q_max,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return ff_calculator
    
    def firin_mah_lazer(self, start_time, end_time, damage_output_handle, target, SPI=False):
        """ 
        end_time: The end time of the photon capture in femtoseconds. Not a real thing experimentally, but useful for choosing 
        a level of damage. Explicitly, it is used to determine the upper time limit for the integration of the form factor.
        pdb_fpath: The pdb file's path. Changes the variable self.atoms.
        y_orientations: 
            Number of unique y axis rotations to sample crystal. x_axis_rotations not implemented (yet?). 
        """
        ff_calculator = self.get_ff_calculator(start_time,end_time,damage_output_handle)     
        target.set_ff_calculator(ff_calculator)
        self.target = target

        ring = np.empty(self.num_rings,dtype="object")
        
        result = Results()
        result.z = 0
        if SPI:
            q_samples = np.linspace(self.q_min,self.q_max,self.num_rings)
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
                        azm = self.alpha_array
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
            # Doing support for single orientation only first. Will need to calculate with new unit cell vectors for each rotation of crystal
            bragg_points = self.bragg_points(target)
            point = np.empty(int(len(bragg_points)),dtype="object")
            radii = np.zeros(len(point))
            azm = np.zeros(len(point))
            z = np.zeros(len(point))
            q_samples = np.zeros(len(point))
            # Get the q vectors where non-zero
            i = 0
            largest_x = -99999
            for G in bragg_points:   # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                if G[0] > largest_x:
                    largest_x = G[0]
                    print("Imaging all G with G[0]",G[0])
                point[i] = self.generate_point(G)
                radii[i] = point[i].r
                azm[i] = point[i].alpha
                q_samples[i] = point[i].q_scr
                z[i] = point[i].I
                #print("G",G,"q",q_samples[i],"z",z[i])
                i+=1
                
            
            r, alph = np.meshgrid(radii, azm) 
            q_for_plot, alph = np.meshgrid(q_samples,azm)  

            #z = np.diag(z)
            # z = np.zeros(r.shape)   # [alpha angle, radius]
            # for i in range(len(point)):
            #     z = point[i].I

            
            
            result.z = z

                # Get the closest angle

                # calculate the 
        
        result.r = r
        result.alph = alph
        result.azm = azm
        result.q_scr = q_for_plot
        self.x_rotation = 0
        return result
    
    class Feature:
        def __init__(self,q,q_scr,r,theta):
            self.q = q
            self.q_scr = q_scr # magnitude of q parallel to screen
            self.r = r
            self.theta = theta          
    class Ring(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    class Spot(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    def generate_ring(self,q):
        '''Returns the intensity(alpha) array and the radius for given q.'''
        #print("q=",q)
        r = self.q_to_r(q)
        q_scr =self.q_to_q_scr(q)
        theta = self.q_to_theta(q)
        ring = self.Ring(q,q_scr,r,theta)
        ring.I = self.illuminate(ring)
        return ring 
    def generate_point(self,G): # G = vector
        q = np.sqrt(G[0]**2+G[1]**2+G[2]**2)
        r = self.q_to_r(q)
        q_scr = self.q_to_q_scr(q)
        theta = self.q_to_theta(q)
        point = self.Spot(q,q_scr,r,theta)
        point.alpha = np.arctan2(G[0],G[1])
        
        alphas = np.array([point.alpha])
        point.I = self.illuminate(point,alphas)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if q_scr != check:
        #     print("Error, q_scr =",q_scr,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity.
    def illuminate(self,feature,alphas = None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet."""
        if alphas == None:
            alphas = self.alpha_array
        F = np.zeros(alphas.shape,dtype="complex_")
        for species in self.target.species_dict.values():
            species.set_scalar_form_factor(feature.q)
            for R in species.coords:
                    # Rotate to crystal's current orientation 
                    R = R.left_multiply(self.y_rot_matrix)  
                    R = R.left_multiply(self.x_rot_matrix)   
                    # Get spatial factor T
                    T = np.zeros(alphas.shape,dtype="complex_")
                    T= self.spatial_factor(alphas,R,feature)
                    F += species.ff*T
                    # Rotate atom for next sample            
        I = np.square(np.abs(F))
        

        return I 

    def spatial_factor(self,alpha_array,R,feature):
        """ theta = scattering angle relative to z-y plane """ 
        q_z = feature.q_scr*np.sin(feature.theta) # not cos because q is hypotenuse - not perpendicular to x-y plane.
        q_z = np.tile(q_z,len(alpha_array))
        # alpha angle of vector relative to x-y plane, pointing from screen centre to point hit. 
        q_y = q_z*np.sin(alpha_array)
        q_x = q_z*np.cos(alpha_array)
        q_vect = np.column_stack([q_x,q_y,q_z])
        for i in range(len(self.target.sym_factors)):
            coord = R.get_array()
            coord = np.multiply(R.get_array(),self.target.sym_factors[i]) + np.multiply(self.target.cell_dim,self.target.sym_translations[i])
            spatial_structure_factor = np.exp(-1j*np.dot(q_vect,coord))   
        return spatial_structure_factor
    def bragg_points(self,crystal, cell_packing = "SC"):
        ''' Using the unit cell structure, find non-zero values of q for which bragg 
        points appear.
        lattice_vectors e.g. = [a,b,c] - length of each spatial vector in orthogonal basis.
        '''
        
        # (h,k,l) is G (subject to selection condition) in lattice vector basis (lattice vector length = 1 in each dimension):
        if cell_packing == "SC":
            lattice_vectors = crystal.cell_dim
        h = 0; k=0; l=0
        b = 2*np.pi/(lattice_vectors)
        # Generate the positive permutations.
        pos_indices = []
        while np.sqrt(((h*b[0])**2+(k*b[1])**2+(l*b[2])**2))<= self.q_max:
            while np.sqrt(((h*b[0])**2+(k*b[1])**2+(l*b[2])**2))<= self.q_max:
                while np.sqrt(((h*b[0])**2+(k*b[1])**2+(l*b[2])**2))<= self.q_max:
                    l0 = l
                    l+=1 
                    if np.sqrt(((h*b[0])**2+(k*b[1])**2)) <= self.q_min:
                        continue                    
                    # Apply Selection conditions. TODO make function.
                    if cell_packing == "SC":  
                        pass
                    if cell_packing == "BCC" and (h+k+l0)%2!=0:
                        continue
                    if cell_packing == "FCC" and (h%2+k%2+l0%2) not in [0,3]:
                        continue
                    pos_indices.append(np.array([h,k,l0]))
                    
                l=0
                k+=1
            k=0
            h+=1     
        h=0       
        # Generate the other permutations(or whatever this is called) (TODO vectorise):
        miller_indices = []
        for h,k,l in pos_indices: 
            sets = [[-h,h],[-k,k],[-l,l]]
            #permutations = np.array(np.meshgrid(h_set,k_set,l_set)).T.reshape(-1, 3)
            #permutations = np.stack(np.meshgrid(h_set,k_set,l_set)).T.reshape(-1, 3)
            permutations = list(itertools.product(*sets))
            permutations = [*set(permutations)]
            miller_indices += permutations
        
        print("Number of points:", len(miller_indices))        
        
        G = np.multiply(b,miller_indices)
        self.rotate_G_to_orientation(G,crystal)
        print(G)
        return G

    def rotate_G_to_orientation(self,G,crystal):
        # TODO not implemented
        return G
                
    # q_abs [1/bohr]
    def q_to_theta(self,q):
        lamb = E_to_lamb(self.photon_energy)
        return np.arcsin(lamb*q/(4*np.pi))   


    def q_to_r(self, q):
        """Returns the distance from centre of screen that is struck by photon which transferred q (magnitude)"""
        D = self.detector_distance #bohr
        theta = self.q_to_theta(q)
        return D*np.tan(2*theta)

    def q_to_q_scr(self,q):
        """ Returns the screen-parallel component of q (or equivalently the final momentum)."""
        theta = self.q_to_theta(q)
        return q*np.cos(2*theta)
    
    def r_to_q(self,r):
        D = self.detector_distance
        lamb = E_to_lamb(self.photon_energy) 
        theta = np.arctan(r/D)/2
        q = 4*np.pi/lamb*np.sin(theta)
        return q 
    def r_to_q_scr(self,r):
        q = self.r_to_q(r)
        return self.q_to_q_scr(q)

    # def q_to_r(self, q_abs):
    #     """Returns the distance from centre of screen that is struck by photon which transferred q"""
    #     D = self.detector_distance #bohr
    #     q0 = self.photon_energy/eV_per_Ha/c_au # initial photon momentum = E/c, directed towards centre of screen.   
    #     theta = self.q_to_theta(q_abs)
    #     qfin_perp = (q0-q_abs)*np.cos(2*theta) # component of the FINAL momentum along beam axis.
    #     r = (D/qfin_perp)*q_abs/np.cos(2*theta)    # -> 0 for q -> 0, -> infinity for 2*theta -> pi/2    # D/qfin_perp is just similar triangle factor.
    #     return r

    # def q_to_q_scr(self,q_abs):
    #     """ Returns the component of the FINAL momentum parallel to the SCREEN."""
    #     theta = self.q_to_theta(q_abs)
    #     q0 = self.photon_energy/eV_per_Ha/c_au
    #     qfin = q0 - q_abs     
    #     print("qfin",qfin)
    #     print("theta",theta)
    #     q_scr = qfin*np.sin(2*theta)
    #     return q_scr
    
    # def r_to_q_scr(self,r):
    #     D = self.detector_distance 
    #     theta = np.arctan(r/D)/2
    #     lamb = E_to_lamb(self.photon_energy)
    #     q_abs = 4*np.pi*np.sin(theta)/lamb
    #     hyp = np.sqrt(D**2+r**2)     
    #     q_scr = (q_abs/hyp)*D               # hyp/D = q/q_scr
    #     # trig check
    #     if q_scr != self.q_to_q_scr(q_abs):
    #         print("error, q_scr:",q_scr,"didn't match what was expected,",self.q_to_q_scr(q_abs)) 
    #     return q_scr  

def E_to_lamb(photon_energy):
    """Energy to wavelength in A.U."""
    E = photon_energy  # eV
    return 2*np.pi*c_au/(E/eV_per_Ha)
    
    # Crystalline.

        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.


def plot_pattern(result,radial_lim = None, plot_against_q=False,log_I = True, log_radial=False,**cmesh_kwargs):
    # https://stackoverflow.com/questions/36513312/polar-heatmaps-in-python
    kw = cmesh_kwargs
    # if "color" not in cmesh_kwargs: 
    #     kw["color"] = 'k'
    # if "ls" not in cmesh_kwargs.keys():
    #     kw["ls"] = 'none'
    fig = plt.figure()
    ax = Axes3D(fig)   

    
    if log_I: 
        z = np.log(result.z)
    else:
        z = result.z
    
    radial_axis = result.r
    if plot_against_q:
        radial_axis =  result.q_scr
    ax = fig.add_subplot(projection="polar")
    if len(result.z.shape) == 1:
        #Point-like
        colours = z
        ax.scatter(result.azm,radial_axis[0],c=colours,**kw)
    else:
        if log_I:
            cutoff_log_intensity = -1
            z -= cutoff_log_intensity
            z[z<0] = 0

        ax.pcolormesh(result.alph, radial_axis, z,**kw)
        ax.plot(result.azm, radial_axis, color = 'k',ls='none')
        plt.grid()  # Make the grid lines represent one unit cell (when implemented).
             
    if radial_lim:
        bottom,top = plt.ylim()
        plt.ylim(bottom,radial_lim)
    if log_radial:
        plt.yscale("log")        
# queue = 0.4
# test = XFEL(12000,100)
# pl_t = test.get_ff_calculator(-10,-9.95,"Naive_Lys_C_7")
# test.set_atomic_species(pl_t,"/home/speno/AC4DC/scripts/scattering/4et8.pdb",["N_fast","S_fast"],CNO_to_N=True)
# print(test.generate_ring(queue).I[0])

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
#%% 
# CNO only.

#energy = 6000 # Tetrapeptide 
energy =17445   #crambin  #q_min=0.11,q_max = 3.9,pixels_per_ring = 400, num_rings = 200
experiment = XFEL(energy,100,x_orientations = 1, y_orientations=1,q_min=0.0175,q_max=3.9, pixels_per_ring = 400, num_rings = 200,t_fineness=100)



sym_translations = [np.array([0,0,0])]
cell_dim = [np.array([1,1,1])]  

experiment.x_rotation = np.pi/2
pdb_path = "/home/speno/AC4DC/scripts/scattering/3u7t.pdb" #Crambin
#pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb"



# 1
allowed_atoms_1 = ["C_fast","N_fast","O_fast"]
end_time_1 = -5
output_handle = "C_tetrapeptide_2"

crystal = Crystal(pdb_path,allowed_atoms_1,CNO_to_N=False)
crystal.set_cell_dim(22.795, 18.826, 41.042)
#crystal.set_cell_dim(10, 10, 10)
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0]))

SPI = False

result1 = experiment.firin_mah_lazer(-10,end_time_1,output_handle,crystal, SPI=SPI)

#%%
# stylin' 
from copy import deepcopy
use_q = False
log_radial = False
log_I = True
cmap = 'Greys'#'binary'
screen_radius = 500
zoom_to_fit = True
####### n'
# plottin'

if zoom_to_fit:
    radial_lim = min(experiment.q_to_r(experiment.q_max),screen_radius)
    print(radial_lim)
if use_q:
    radial_lim = experiment.q_to_q_scr(experiment.q_max)
    print(radial_lim)

result_mod = deepcopy(result1)
result_mod.z = result_mod.z #hack to remove neg nums (due to taking log).

plot_pattern(result_mod,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I)
#%% 2
allowed_atoms_2 = ["C_fast","N_fast","O_fast"]
end_time_2 = -5
output_handle = "C_tetrapeptide_2"
result2 = experiment.firin_mah_lazer(-10,end_time_2,output_handle,pdb_path,allowed_atoms_2,CNO_to_N=False)
experiment.plot_pattern(result2,plot_against_q = use_q,log=use_log)
#%% Difference
result3 = Results()
result3.r = result1.r
result3.q = result1.q
result3.z = np.abs(np.sqrt(np.exp(result1.z))-np.sqrt(np.exp(result2.z)))/np.sqrt(np.exp(result1.z))
print(result3.z[0][0])
print(np.abs(np.sqrt(np.exp(result1.z[0][0]))-np.sqrt(np.exp(result2.z[0][0])))/np.sqrt(np.exp(result1.z[0][0])))
result3.alph = result1.alph
result3.azm = result1.azm
experiment.plot_pattern(result3,plot_against_q = use_q,log=use_log)

#%% Tetra experimental conditions kinda 
detector_dist = 100
experiment = XFEL(12562,detector_dist,x_orientations = 1, y_orientations=1, pixels_per_ring = 200, num_rings = 250,t_fineness=100)
#experiment.q_max = experiment.r_to_q(detector_dist*1.99) 
experiment.q_max = 1/(1.27*1.8897) 
experiment.q_min = 1/(14.83*1.8897)  
use_q = False
use_log = False
pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb"

# 1

allowed_atoms_1 = ["C_fast"]
end_time_1 = -9.99
dmg_output = "C_tetrapeptide_2"
result1 = experiment.firin_mah_lazer(-10,end_time_1,dmg_output,pdb_path,allowed_atoms_1,CNO_to_N=False)
#%%

experiment.plot_pattern(result1,plot_against_q = use_q,log=use_log)
plt.yscale("symlog")

#TODO 
# log doesnt work atm.
# %%


# get sinc functions with non-infinite crystal
# 