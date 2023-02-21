#%%
#TODO make class
#%%
import os
os.getcwd() 
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
print(sys.path)

from Bio.PDB.vectors import Vector as bio_vect
from Bio.PDB.PDBParser import PDBParser
import numpy as np

import matplotlib.pyplot as plt
from plotter_core import Plotter

# Get what we need
handle = "Naive_Lys_C_7"
parser=PDBParser(PERMISSIVE=1)
structure_id="4et8"
filename="/home/speno/AC4DC/scripts/pdb_parser/4et8.pdb"

pl = Plotter(handle,"y")
plt.close()
photon_energy = 6000
q_max = 2 # In units of bohr^-1. 
pl.initialise_coherence_params(-10,-9.6,q_max,photon_energy,50,100,True)


#%%
do_print = True

def spatial_factor_but_wrong(q,R):
    global do_print
    # Spherical averaging
    sq_sample_size = 100  # 
    # Scattering angle
    theta = np.linspace(0,np.pi,sq_sample_size,endpoint=False)
    # Angle of vector relative to arbitrary angle, pointing from screen centre to point hit. 
    alpha = np.linspace(0,2*np.pi,sq_sample_size,endpoint=False)
    
    q_z = q*np.sin(theta) # not cos because q is hypotenuse - not perpendicular to x-y plane.
    q_z = np.tile(q_z,(sq_sample_size,1)) # copies of rows containing z factor for each sample. Number of copies corresponds to number of alpha angles we sample
    y_alph = np.tile(np.sin(alpha),(sq_sample_size,1)).transpose()
    q_y = np.multiply(q_z,y_alph) 
    x_alph = np.tile(np.cos(alpha),(sq_sample_size,1)).transpose()
    q_x = np.multiply(q_z,x_alph)


    q_x = q_x.flatten()
    q_y = q_y.flatten()
    q_z = q_z.flatten()
    q_vect = np.column_stack([q_x,q_y,q_z])  

    R = R.get_array() * 1.88973  # convert to bohr
    
    neg_i_q_R = -1j*np.dot(q_vect,R)
    samples = np.exp(neg_i_q_R)

    
    return np.average(samples)
    
    ############ non-vectorised
    '''
    # Spherical averaging
    sq_sample_size = 4
    # Scattering angle
    theta = np.linspace(0,np.pi,sq_sample_size,endpoint=False)
    # Angle of vector relative to arbitrary angle, pointing from screen centre to point hit. 
    alpha = np.linspace(0,2*np.pi,sq_sample_size,endpoint=False)
    
    samples = np.zeros(sq_sample_size**2,dtype="complex_")
    for j,t in enumerate(theta):
        q_z = q*np.sin(t) # not cos because q is hypotenuse - not perpendicular to x-y plane.
        for i,a in enumerate(alpha):
            q_y = q_z*np.sin(a)
            q_x = q_z*np.cos(a)
            q_vect = bio_vect(q_x,q_y,q_z)
            samples[i*sq_sample_size+j] = np.exp(-1j*np.dot(q_vect,R))    
    print(np.average(samples))
    return np.average(samples)    
    '''


# focused beam --> theta we know, but have range of alpha.
def spatial_factor(q, theta, R):
    # theta = scattering angle relative to z-y plane
    z_comp = q*np.sin(theta) # not cos because q is hypotenuse - not perpendicular to x-y plane.
    sample_size = 400
    alpha = np.linspace(0,2*np.pi,sample_size,endpoint=False)
    # alpha angle of vector relative to x-y plane, pointing from screen centre to point hit. 
    q_y = z_comp*np.sin(alpha)
    q_x = z_comp*np.cos(alpha)
    q_z = np.zeros(sample_size)
    q_z.fill(z_comp)
    q_vect = np.column_stack([q_x,q_y,q_z])
    samples = np.exp(-1j*np.dot(q_vect,R.get_array()))    
    return np.average(samples)   

# Polarised light
def spatial_factor_polarised(q, theta, R, alpha):
    # theta = scattering angle relative to z-y plane
    q_z = q*np.sin(theta) # not cos because q is hypotenuse - not perpendicular to x-y plane.
    # alpha angle of vector relative to x-y plane, pointing from screen centre to point hit. 
    q_y = q_z*np.sin(alpha)
    q_x = q_z*np.cos(alpha)
    q_vect = np.column_stack([q_x,q_y,q_z])
    return np.exp(-1j*np.dot(q_vect,R.get_array()))   

def illuminate_polarised(q, theta, atoms, alpha_array,simple = True):
    def fast_name(name):
        allowed_atoms = ["N","S"]
        if simple:
            if name == "C" or name == "O":
                name = "N"
            if name not in allowed_atoms:
                return None
        if name != "H":
            name+="_fast"
        return name
    F = np.zeros(alpha_array.shape,dtype="complex_")
    for atom in atoms:
        name = atom.element
        a = fast_name(name)
        if a == None:
            continue
        # F_i(q) = f(q)*T_i(q), where f is the time-integrated average 
        f = pl.f_average(q,a)
        R = atom.get_vector()
        T = np.zeros(alpha_array.shape,dtype="complex_")
        for a,alpha in enumerate(alpha_array):  #TODO vectorise
            T[a]= spatial_factor_polarised(q,theta,R,alpha)
        F += f*T
    I = np.square(np.abs(F))

    return I


# Crystalline.
def illuminate_miller(atom):
    pass
    #non-zero q: miller indices q_hkl

# Get the contribution to the intensity term  from the atom. Not crystalline.
def illuminate(q, theta, atoms,simple = True):
    def fast_name(name):
        allowed_atoms = ["N","S"]
        if simple:
            if name == "C" or name == "O":
                name = "N"
            if name not in allowed_atoms:
                return None
        if name != "H":
            name+="_fast"
        return name
    F = 0
    for atom in atoms:
        name = atom.element
        a = fast_name(name)
        if a == None:
            continue
        # F_i(q) = f(q)*T_i(q), where f is the time-integrated average 
        f = pl.f_average(q,a)
        R = atom.get_vector()
        T= spatial_factor(q,theta,R)
        F += f*T
    I = np.abs(F)**2

    return I

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

# Returns the intensity at q
def intensity(q,theta,alpha_array):
    atoms = []
    structure=parser.get_structure(structure_id, filename)
    model=structure[0]
    for chain in model.get_list():
        atoms_gen = chain.get_atoms()
        for atom in atoms_gen:
            atoms.append(atom)
    if alpha_array is not None:
        return illuminate_polarised(q,theta,atoms,alpha_array)
    return illuminate(q, theta, atoms)


# E = photon energy in eV
def q_to_theta(q_abs,E):
    bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h_ev = 6.582119569e-16 *2 *np.pi  #TODO check why angstrom_per_nm is used instead by HR but with bohr units.
    lamb =  h_ev*c_nm/E
    lamb *= bohr_per_nm   
    return np.arcsin(lamb*q_abs/(4*np.pi))    

def q_to_x(q_abs,E,D):
    theta = q_to_theta(q_abs,E)
    return D*np.tan(theta)    

# def x_to_q(x,E,D):
#     # D*np.tan(np.arcsin(lamb*q_abs/(4*np.pi))) = x, q = sin(arctan((x/D)))/lamb*4*np.pi
#     bohr_per_nm = 18.8973; c_nm = 2.99792458e17; h_ev = 6.582119569e-16 *2 *np.pi  #TODO check why angstrom_per_nm is used instead by HR but with bohr units.
#     lamb =  h_ev*c_nm/E
#     lamb *= bohr_per_nm       
#     q = np.sin(np.arctan((x/D)))/lamb*4*np.pi
#     return q

# Returns the intensity at a position on the screen and the x value for given q.
def scattering_pattern(q,E,D,alpha_array=None):
    x = q_to_x(q,E,D)
    theta = q_to_theta(q,E)
    print("q=",q)
    I = intensity(q,theta,alpha_array) 
    return x,I

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
print(queue,scattering_pattern(queue,12000,100))


#%%
handle = "Improved_Lys_early_11"
parser=PDBParser(PERMISSIVE=1)
structure_id="4et8"
filename="/home/speno/AC4DC/scripts/pdb_parser/4et8.pdb"

#%% Polarised EARLIER
E_TIME = -9.95
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from plotter_core import Plotter

pl = Plotter(handle,"y")
plt.close()
photon_energy = 6000
q_max = 2.4 # In units of bohr^-1. 
pl.initialise_coherence_params(-10,E_TIME,q_max,photon_energy,50,100,True)

q_size = 100
alpha_size = 400
y = np.zeros(q_size,dtype="object")
R = np.zeros(q_size,dtype="float")
alpha_array = np.linspace(0,2*np.pi,alpha_size,endpoint=False)
q_samples = np.linspace(0.5,2.4,q_size)
for i, q in enumerate(q_samples):
    R[i], y[i] = scattering_pattern(q,12000,100,alpha_array)
    y[i] = (y[i]+1)
    print("q:",q, "x:",R[i],"I[alph=0]",y[i][0])

fig = plt.figure()
ax = Axes3D(fig)

azm = alpha_array
r, alph = np.meshgrid(R, azm)
z = np.zeros(r.shape)  # z is the intensity of the plot colour.
for pos in range(len(z)):
    for ang in range(len(z[pos])):
        z[pos][ang] = np.log(y[ang][pos]) 

print(r)
print(z)

plt.subplot(projection="polar")
plt.pcolormesh(alph, r, z)
plt.plot(azm, r , color='k', ls='none') 
plt.grid()  # Make the grid lines represent one unit cell when we do that. 
plt.show()

#%% LATER
L_TIME = -9.7
pl = Plotter(handle,"y")
plt.close()
q_max = 2.4 # In units of bohr^-1. 
pl.initialise_coherence_params(-10,L_TIME,q_max,photon_energy,50,100,True)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
q_size = 100
alpha_size = 400
y2 = np.zeros(q_size,dtype="object")
R2 = np.zeros(q_size,dtype="float")
alpha_array2 = np.linspace(0,2*np.pi,alpha_size,endpoint=False)
q_samples2 = np.linspace(0.5,q_max,q_size) 
for i, q in enumerate(q_samples2):
    R2[i], y2[i] = scattering_pattern(q,12000,100,alpha_array2)
    y2[i] = (y2[i]+1)
    print("q:",q, "x:",R2[i],"I[alph=0]",y2[i][0])

fig2 = plt.figure()
ax2 = Axes3D(fig2)

azm2 = alpha_array2
r2, alph2 = np.meshgrid(R2, azm2)
z2 = np.zeros(r2.shape)  # z is the intensity of the plot colour.
for pos in range(len(z2)):
    for ang in range(len(z2[pos])):
        z2[pos][ang] = np.log(y2[ang][pos]) 

print(r2)
print(z2)

plt.subplot(projection="polar")
plt.pcolormesh(alph2, r2, z2)
plt.plot(azm2, r2 , color='k', ls='none') 
plt.grid()
plt.show()



#%% COMPARE
z_diff = z-z2
alph3 = alph
r3 = r

print(z_diff)
fig3 = plt.figure()
ax3 = Axes3D(fig2)

plt.subplot(projection="polar")
plt.pcolormesh(alph3, r3, z_diff)
plt.plot(alph3, r3 , color='k', ls='none') 



# %% Not polarised
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
size = 50
y = np.zeros(size)
R = np.zeros(size)
alpha = 0
q_samples = np.linspace(0.5,2.4,size)
for i, q in enumerate(q_samples):
    R[i], y[i] = scattering_pattern(q,12000,100)
    y[i] = (y[i]+1)
    print("q:",q, "x:",R[i],"I",y[i])

fig = plt.figure()
ax = Axes3D(fig)

azm = np.linspace(0, 2*np.pi, size)
r, th = np.meshgrid(R, azm)
z = (r ** 2.0) / 4.0
for i in range(len(z)):
    for j in range(len(z[i])):
        z[i][j] = np.log(y[j])

print(r)
print(z)

plt.subplot(projection="polar")
plt.pcolormesh(th, r, z)
plt.plot(azm, r , color='k', ls='none') 
plt.grid()
plt.show()
# %%


fig = plt.figure()
ax = Axes3D(fig)

rad = x
azm = np.linspace(0, 2 * np.pi, 100)
r, th = np.meshgrid(rad, azm)
z = (r ** 2.0) / 4.0

plt.subplot(projection="polar")

plt.pcolormesh(th, r, z)
#plt.pcolormesh(th, z, r)

plt.plot(azm, r, color='k', ls='none') 
plt.grid()

plt.show()
# %%
