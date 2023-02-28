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
    def __init__(self, pdb_fpath, allowed_atoms, cell_packing = "SC", CNO_to_N = False, orbitals_as_shells = True):
        self.cell_packing = cell_packing
        
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
    def __init__(self, photon_energy, detector_distance,  x_orientations = 1, y_orientations = 1,q_max=2.4, pixels_per_ring = 400, num_rings = 50,t_fineness=100):
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
            q_samples = np.linspace(0,self.q_max,self.num_rings)
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
            bragg_points, miller_indices = self.bragg_points(target,cell_packing =target.cell_packing)
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
        if not SPI:
            result.txt = miller_indices
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
        print(" ")
        print(np.sqrt(G[0]**2+G[1]**2))
        print(self.q_to_q_scr(q))
        r = self.q_to_r(q)
        q_scr = np.sqrt(G[0]**2+G[1]**2)#self.q_to_q_scr(q)#np.sqrt(G[0]**2+G[1]**2) #self.q_to_q_scr(q)  # not np.sqrt(G[0]**2+G[1]**2)
        theta = self.q_to_theta(q)
        point = self.Spot(q,q_scr,r,theta)
        point.alpha = np.arctan2(G[1],G[0])
        point.G = G
        print(point.alpha)

                

        print("q_scr",q_scr,"alpha",point.alpha,"G",G[0],G[1],G[2])


        #print("alpha",point.alpha,"q_scr",point.q_scr)
        
        alphas = np.array([point.alpha])
        point.I = self.illuminate(point,alphas)
        # # Trig check      
        # check = smth  # check = np.sqrt(G[0]**2+G[1]**2)
        # if q_scr != check:
        #     print("Error, q_scr =",q_scr,"but expected",check)  
        return point   
    

    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
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
    def bragg_points(self,crystal, cell_packing):
        ''' Using the unit cell structure, find non-zero values of q for which bragg 
        points appear.
        lattice_vectors e.g. = [a,b,c] - length of each spatial vector in orthogonal basis.
        '''
        
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
        q_x_max = self.q_max#self.q_to_q_scr(self.q_max)   # Max length of any vector parallel to screen.
        q_y_max = q_x_max
        q_z_max = self.q_max        # q = (0,0,l) case.
        h_max = 0
        k_max = 0
        l_max = 0
        h = 0
        while h*np.sqrt(sum(pow(element, 2) for element in b[0])) <= q_x_max:
            h_max += 1
            h += 1
        k = 0
        while k*np.sqrt(sum(pow(element, 2) for element in b[1])) <= q_y_max:
            k_max += 1
            k += 1      
        l = 0
        while l*np.sqrt(sum(pow(element, 2) for element in b[2])) <= q_z_max:
            l_max += 1
            l += 1                     
        h = 0; k = 0;l = 0

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

        G_temp = np.zeros(indices.shape)
        for i in range(len(indices)):
            G_temp[i] = np.array(np.dot(indices[i],b))
        self.rotate_G_to_orientation(G_temp,crystal) 

        # Purely a matter of truncation
        q_max_rule = lambda g: np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2))<= self.q_max
        #q_max_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.q_max
        # Since our screen is perpendicular to the incoming beam, we have an additional constraint from elasticity:       
        elasticity_constraint = lambda g: ((g[2] > 0 and (g[0] != 0 or g[1] != 0)) or (g[2] == 0 and g[0] == 0 and g[1] == 0))
        
        # Catch the miller indices *cough cou-gh-ish* with a boolean mask
        mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(q_max_rule,1,G_temp)*np.apply_along_axis(elasticity_constraint,1,G_temp)
        indices = indices[mask]

        print("Number of points:", len(indices))   
        for elem in indices:
            print(elem)
        
        #TODO vectorise somehow
        G = np.zeros(indices.shape)
        for i in range(len(indices)):
            G[i] = np.array(np.dot(indices[i],b))
        self.rotate_G_to_orientation(G,crystal)
        print(G)
        return G, indices

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
        """ Returns the screen-parallel component of q (or equivalently the final momentum).
        This whole function is dubious. Need to replace probably.
        """
        theta = self.q_to_theta(q)
        return q*np.cos(theta)
    
    def r_to_q(self,r):
        D = self.detector_distance
        lamb = E_to_lamb(self.photon_energy) 
        theta = np.arctan(r/D)/2
        q = 4*np.pi/lamb*np.sin(theta)
        return q 
    def r_to_q_scr(self,r):
        q = self.r_to_q(r)
        return self.q_to_q_scr(q)

def E_to_lamb(photon_energy):
    """Energy to wavelength in A.U."""
    E = photon_energy  # eV
    return 2*np.pi*c_au/(E/eV_per_Ha)
    
    # Crystalline.

        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.


def scatter_plot(result,show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None, **cmesh_kwargs):
    '''When crystalline, dot size is proprtional to intensity, while colour is proportional to natural log of intensity.'''
    # https://stackoverflow.com/questions/36513312/polar-heatmaps-in-python
    kw = cmesh_kwargs
    # if "color" not in cmesh_kwargs: 
    #     kw["color"] = 'k'
    # if "ls" not in cmesh_kwargs.keys():
    #     kw["ls"] = 'none'
    fig = plt.figure()
    ax = Axes3D(fig)   

       
    radial_axis = result.r
    if plot_against_q:
        radial_axis =  result.q_scr
    ax = fig.add_subplot(projection="polar")
    if len(result.z.shape) == 1:
        radial_axis = radial_axis[0]
        ## Point-like (Crystalline)
        # get points at same position on screen.
        identical_count = np.zeros(result.z.shape)
        tmp_z = np.zeros(result.z.shape)
        processed_copies = []
        unique_values_mask = np.zeros(result.z.shape,dtype="bool")
        for i in range(len(result.z)):
            if (radial_axis[i], result.azm[i])  in processed_copies:
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
            processed_copies.append((radial_axis[i],result.azm[i]))
        
        identical_count = identical_count[unique_values_mask]
        tmp_z = tmp_z[unique_values_mask]
        radial_axis = radial_axis[unique_values_mask]
        azm = result.azm[unique_values_mask]
        
        z = tmp_z # Idk why I had this: #np.multiply(tmp_z,identical_count)

        if log_I: 
            z = np.log(z)
        else:
            z = z  

        debug_mask = (0.01 < radial_axis[0])*(radial_axis[0] < 100)
        print(identical_count[debug_mask])
        print(radial_axis[0][debug_mask])
        print(azm[debug_mask ]*180/np.pi)

        colours = z
        # Dot size
        dot_param = z
        if crystal_pattern_only:
            dot_param = identical_count
        if not log_dot:
            dot_param = np.e**(dot_param)     
        norm = np.max(dot_param)
        s = [100*dot_size*x/norm for x in dot_param]

        ax.scatter(azm,radial_axis,c=colours,s=s,alpha=1,**kw)
        plt.grid() 
        if show_labels:
            for i, txt in enumerate(result.txt[unique_values_mask]):
                ax.annotate("%.0f" % txt[0]+","+"%.0f" % txt[1]+","+"%.0f" % txt[2]+",", (azm[i], radial_axis[i]))        
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
experiment.scatter_plot(result1)
#%%
allowed_atoms_2 = ["N_fast","S_fast"]
end_time_2 = -9.76
output_handle = "Naive_Lys_C_7"
pdb_path = "/home/speno/AC4DC/scripts/scattering/4et8.pdb"
result2 = experiment.firin_mah_lazer(-10,end_time_2,output_handle,pdb_path,allowed_atoms_2,CNO_to_N=True)
experiment.scatter_plot(result2)

#%% Difference
result3 = Results()
result3.r = result1.r
result3.q = result1.q
result3.z = result1.z-result2.z
result3.alph = result1.alph
result3.azm = result1.azm
experiment.scatter_plot(result3)
#
#%% 
# CNO only.

#energy = 6000 # Tetrapeptide 
energy =17445   #crambin - from resolutions. in pdb file, need to double check calcs: #q_min=0.11,q_max = 3.9,  my defaults: pixels_per_ring = 400, num_rings = 200
experiment = XFEL(energy,100,x_orientations = 1, y_orientations=1,q_max=9, pixels_per_ring = 400, num_rings = 400,t_fineness=100)

experiment.x_rotation = 0#np.pi/2#0#np.pi/2

sym_translations = [np.array([0,0,0])]
cell_dim = [np.array([1,1,1])]  

pdb_path = "/home/speno/AC4DC/scripts/scattering/3u7t.pdb" #Crambin
#pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb"



# 1
#allowed_atoms_1 = ["C_fast","N_fast","O_fast"]
allowed_atoms_1 = ["C_fast","N_fast","O_fast","S_fast"]
#end_time_1 = -5
#output_handle = "C_tetrapeptide_2"
end_time_1 = -9.8
output_handle = "Improved_Lys_mid_6"

crystal = Crystal(pdb_path,allowed_atoms_1,CNO_to_N=False,cell_packing = "BCC")
#crystal.set_cell_dim(22.795*1.88973, 18.826*1.88973, 41.042*1.88973)
crystal.set_cell_dim(20, 20, 20)
crystal.add_symmetry(np.array([-1, 1,-1]),np.array([0,0.5,0]))

SPI = False

result1 = experiment.firin_mah_lazer(-10,end_time_1,output_handle,crystal, SPI=SPI)

#%%
# stylin' 

font = {'family' : 'monospace',
        'weight' : 'bold',
        'size'   : 24}

plt.rc('font', **font)

from copy import deepcopy
use_q = True
log_radial = False
log_I = True
cutoff_log_intensity = -1#-1
cmap = 'Greys'#'binary'
#screen_radius = 1400
screen_radius = 120#55#165    #
q_scr_lim = experiment.r_to_q_scr(screen_radius)#3.9
zoom_to_fit = True
####### n'
# plottin'

if zoom_to_fit:
    radial_lim = min(screen_radius,experiment.q_to_r(experiment.q_max))
    print(radial_lim)
    if use_q:
        #radial_lim = experiment.r_to_q_scr(radial_lim)
        radial_lim = min(q_scr_lim,experiment.q_to_q_scr(experiment.q_max))
        print(radial_lim)
else:
    radial_lim = None

result_mod = deepcopy(result1)
print(result_mod.z)
result_mod.z = result1.z
print(np.log(result_mod.z))

radial_lim = 10
scatter_plot(result_mod,crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity)

fig = plt.gcf()
fig.set_figwidth(20)
fig.set_figheight(20)

#NEED TO CHECK. We have a 1:1 mapping from q to q_parr, but with our miller indices we are generating multiple q_parr with diff q.
# So we SHOULD get the same q_parr with different q_z. Which makes sense since we are just doing cosine. But still icky maybe?
# Need to double check we get different intensities for same q_parr. Pretty sure that's implemented.
#%% DEBUG 
#print(np.arctan2(3,-1))
#print(experiment.q_to_theta(np.sqrt(11)))
#print(experiment.q_to_q_scr(np.sqrt(11)))
print(result_mod.z)
print(result_mod.r[0])
print(result_mod.q_scr[0])
#%% 2
allowed_atoms_2 = ["C_fast","N_fast","O_fast"]
end_time_2 = -5
output_handle = "C_tetrapeptide_2"
result2 = experiment.firin_mah_lazer(-10,end_time_2,output_handle,pdb_path,allowed_atoms_2,CNO_to_N=False)
experiment.scatter_plot(result2,plot_against_q = use_q,log=use_log)
#%% Difference
result3 = Results()
result3.r = result1.r
result3.q = result1.q
result3.z = np.abs(np.sqrt(np.exp(result1.z))-np.sqrt(np.exp(result2.z)))/np.sqrt(np.exp(result1.z))
print(result3.z[0][0])
print(np.abs(np.sqrt(np.exp(result1.z[0][0]))-np.sqrt(np.exp(result2.z[0][0])))/np.sqrt(np.exp(result1.z[0][0])))
result3.alph = result1.alph
result3.azm = result1.azm
experiment.scatter_plot(result3,plot_against_q = use_q,log=use_log)

#%% Tetra experimental conditions kinda 
detector_dist = 100
experiment = XFEL(12562,detector_dist,x_orientations = 1, y_orientations=1, pixels_per_ring = 200, num_rings = 250,t_fineness=100)
#experiment.q_max = experiment.r_to_q(detector_dist*1.99) 
experiment.q_max = 1/(1.27*1.8897) 
use_q = False
use_log = False
pdb_path = "/home/speno/AC4DC/scripts/scattering/5zck.pdb"

# 1

allowed_atoms_1 = ["C_fast"]
end_time_1 = -9.99
dmg_output = "C_tetrapeptide_2"
result1 = experiment.firin_mah_lazer(-10,end_time_1,dmg_output,pdb_path,allowed_atoms_1,CNO_to_N=False)
#%%

experiment.scatter_plot(result1,plot_against_q = use_q,log=use_log)
plt.yscale("symlog")

#TODO 
# log doesnt work atm.
# %%


# get sinc functions with non-infinite crystal
# 