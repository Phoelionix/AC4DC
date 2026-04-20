#%%
#TODO
## Important
# - make it so reflections don't overwrite same orientation, as stochastic now.
# - should have option to average out same miller indices be averaged out.
# - normalise damaged I rather than using neutze-style k factor in R factor.
## Not so important 
# - implement rhombic miller indices as the angle is actually 120 degrees on one unit cell lattice vector (or just do SPI)
# - Select times and integrate snapshots of intensities with gaussian quadrature or some better method than trapezoid method. (I started this with scatter_quad.py)
# - Make all paths absolute

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

import os
os.getcwd()
import sys
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/')
######

import os.path as path
from Bio.PDB.vectors import Vector as Bio_Vect
from Bio.PDB.vectors import homog_trans_mtx, set_homog_trans_mtx
from Bio.PDB.vectors import rotaxis2m
#from Bio.PDB.PDBParser import PDBParser
#from Bio.PDB.PDBIO import PDBIO
#from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.MMCIFParser import MMCIFParser # no xpdb verion
from xpdb import sloppyparser as xPDBParser
from xpdb import SloppyPDBIO as xPDBIO
from xpdb import SloppyStructureBuilder as xStructureBuilder  # Hack, enables atom counts over 10,000
from Bio.PDB.Atom import Atom as PDB_Atom
from Bio.PDB.Atom import DisorderedAtom
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
from matplotlib.colors import TwoSlopeNorm
import copy
import pickle
#import colorcet as cc; import cmasher as cmr
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.offline as pltly_offline
from IPython.display import display, HTML
from IPython import get_ipython
from string import ascii_uppercase, ascii_lowercase, ascii_letters, digits
from core_functions import get_sim_elements,ATOMNO,parse_elecs_from_latex
from plot_I_vs_Isigma import read_cif
import scipy
import subprocess
from contextlib import contextmanager,redirect_stderr,redirect_stdout
from os import devnull
import labellines
from multiprocessing import Pool

NUM_THREADS=20

interactive = True
if interactive and __name__ == "__main__":
    get_ipython().run_line_magic('colors', 'nocolor')
    get_ipython().run_line_magic('matplotlib', 'widget')
    # pltly_offline.init_notebook_mode()
    # display(HTML(
    #     '<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_SVG"></script>'
    # ))    

plt.ioff()  # stops weird vscode stuff

USE_PHENIX=True

DEBUG = False; DEBUG = False; DEBUG_MODERATE = False; DEBUG_WATER = True
SEEDED = False# TODO check fully implemented for all random stuff
DELETENONUNIQUE=False

if SEEDED:
    np.random.seed(0)

c_au = 137.036; eV_per_Ha = 27.211385; ang_per_bohr = 1/1.8897259886 # > 1 ang = 1.88973 bohr

RESULTS_LOCAL_PATH = "results/"

@contextmanager
# def suppress_stdout():
#     #https://stackoverflow.com/questions/2125702/how-to-suppress-console-output-in-python
#     with open(os.devnull, "w") as devnull:
#         old_stdout = sys.stdout
#         sys.stdout = devnull
#         try:  
#             yield
#         finally:
#             sys.stdout = old_stdout
def suppress_stdout_stderr():
    #https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions
    with open(devnull, 'w') as fnull:
        with redirect_stderr(fnull) as err, redirect_stdout(fnull) as out:
            #yield (err, out)
            yield (err, out)
    
class Custom_Gromacs_Parser():
    class Structure():
        def __init__(self,structure_id,conf_path):
            self.id = structure_id
            self.conf_path = conf_path
        def get_atoms(self):
            atoms = []
            NA_warning=False; CL_warning=False
            with open(self.conf_path) as gromacs_config_file:
                i = 1
                header_remaining = 2
                for line in gromacs_config_file:
                    # Figure out where the heck we are. May or may not work for all .gro files.
                    #####
                    if header_remaining > 0:
                        header_remaining -= 1
                        continue
                    if "." in line.split()[0]: # generally the last line with cell dimensions
                        continue
                ######
                    assert len(line)>20
                    vals = line[0:8].strip(), line[11:15].strip(), line[15:20].strip(), *(line[20:].split())
                    # After 10k the name and serial number columns are joined together
                    #elif len(vals[1]) > 4: # NB: If reach  99999, then the index resets.
                        # vals[2] = last_val+1
                        # oom = len(str(i%1e5)) # NB: If reach  99999, then the index resets.
                        # vals.insert(2,vals[1][-oom:])
                        # vals[1] = vals[1][:-oom]


                    element = None    
                        # NOTE if modify here, need to modify in parser
                    name = vals[1].strip()

                    # Commented out to avoid confusing with atom labelled NA in HEME.
                    # TODO put in pdb warning
                    if name == "NA":
                        NA_warning = True
                    if name == "CL":
                        CL_warning=True
                        #element = name
                    special_convert_dict={"FE":"FE","CLA":"CL","SOD":"NA","ZN":"ZN","CAL":"CA"} # changes here should be made below
                    if name in special_convert_dict:
                        element = special_convert_dict[name]
                    with suppress_stdout_stderr():
                        atom = PDB_Atom(
                            name = vals[1],
                            coord = (float(vals[3])*10,float(vals[4])*10,float(vals[5])*10), # converts from nm to angstrom
                            bfactor = 0,
                            occupancy = None,
                            altloc = None,
                            fullname = " " + vals[1] + " ",
                            serial_number = i,
                            element=element
                        )
                    atoms.append(atom)
                    assert(int(vals[2])==i%1e5), (vals[2],i,"||", line, "||", vals,"||",last_val,last_i)
                    last_val, last_i = int(vals[2]),i
                    i+=1
            if NA_warning:
                print(f"Warning: atom with name NA will be treated as nitrogen")
            if CL_warning:
                print(f"Warning: atom with name CL will be treated as carbon")
            return atoms

    def get_structure(self,structure_id, conf_path):
            return self.Structure(structure_id,conf_path)

class Results():
    def __init__(self,target:'Crystal',num_points=None,image_index=None):
        self.cell_dims= [v*ang_per_bohr for v in target.cell_dim]
        self.cell_angles=target.cell_angles
        self.symmetry=target.symmetry
        
        self.phi = np.zeros(num_points)
        self.phi_aligned = np.zeros(num_points)
        self.I = np.zeros(num_points)
        self.q = np.zeros(num_points)
        self.X = np.zeros(num_points)
        self.image_index = image_index 
        self.for_plotting=True

        
    @staticmethod
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

    def package_up(self,miller_indices,for_plotting=True):
        # _, self.phi_mesh = np.meshgrid(self.q,self.phi)  
        # self.X, self.phi_mesh = np.meshgrid(self.X, self.phi)
        # self.q, self.phi_aligned_mesh = np.meshgrid(self.q,self.phi_aligned) 
        self.for_plotting=for_plotting
        if for_plotting:
            self.X, _ = np.meshgrid(self.X, self.phi)
            self.q, _ = np.meshgrid(self.q,self.phi_aligned) 
        self.miller_indices = miller_indices 
    def diff(self,other):
        '''
        Get single-point R factors between two images
        '''
        self.R = np.zeros(self.I.shape)

        if not self.for_plotting:
            return

        
        CHECK_ALIGNED = True
        if CHECK_ALIGNED:
            for i in range(len(self.q)):
                subtracted = False
                for j in range(len(other.q)):
                    if self.phi[i] == other.phi[j] and self.q[0][i] == other.q[0][j]:   # Using phi is equivalent to using phi_aligned, as this function is only used for same-orientation comparisons.
                        #self.I[i] = np.abs(self.I[i]- other.I[j])/self.I[i] # Intensity
                        self.R[i] = np.abs(np.sqrt(self.I[i])- np.sqrt(other.I[j]))/np.sqrt(max(self.I[i],other.I[j]))# form factor  # TODO not sure why later times is sometimes larger but probs interference thing. Hopefully will disappear when we use large time gaps... or at least when we average over to get R factor.
                        subtracted = True
                if subtracted == False:
                    self.R[i] = -1 
        else:
            for i in range(len(self.q)):
                self.R[i] = np.abs(np.sqrt(self.I[i])- np.sqrt(other.I[i]))/np.sqrt(max(self.I[i],other.I[i]))

class Results_SPI():
    def package_up(self,for_plotting=True):
        self.for_plotting=for_plotting

class Crystal():
    def __init__(self, struct_file_path, allowed_atoms, positional_stdv = 0, is_damaged=True, include_symmetries = None, rocking_angle = 0.3, 
    cell_packing = "SC", CNO_to_N = False, supercell_scale = 1,num_supercells=1, supercell_simulations = 1, 
    S_to_N=False,convert_excluded_elements_to_H=False,convert_excluded_elements_to_N=False,allow_skip_species=False,random_waters=None,
    use_bfactors=True,zero_bfactors=False,ignore_water_H=False,charge_states=None, use_intensity_for_time=None):
        allowed_atoms=copy.deepcopy(allowed_atoms)
        '''
        rocking_angle [degrees]
        cell_packing ("SC","BCC","FCC","FCC-D", "triclinic")
        '''
        self.gromacs_config_file = struct_file_path.split('.')[-1]=="gro"
        self.cif_file = struct_file_path.split('.')[-1]=="cif"
        assert not self.cif_file, "cif not supported"
        if include_symmetries is None:
            include_symmetries = not self.gromacs_config_file
        if self.gromacs_config_file:
            print("Using gromacs file")
            assert include_symmetries == False
            assert num_supercells == 1
            
        if zero_bfactors:
            assert use_bfactors, "Can't set zero B factors - B factors aren't being used."

        if use_intensity_for_time is not None:
            assert not is_damaged


        self.stochastic_positions_set=False
        
        assert not (convert_excluded_elements_to_H and convert_excluded_elements_to_N)

        self.use_bfactors = use_bfactors
        self.zero_bfactors = zero_bfactors
        self.cell_packing = cell_packing
        self.rocking_angle = rocking_angle * np.pi/180            
        self.is_damaged = is_damaged
        self.supercell_scale = supercell_scale  # TODO allow for non-cubic crystals and non SC cell packing.
        self.num_supercells = int(num_supercells)
        self.supercell_simulations = supercell_simulations
        if positional_stdv == 0 and is_damaged == False:
            self.supercell_simulations = 1
        self.struct_file_path = struct_file_path
        self.positional_stdv = positional_stdv/ang_per_bohr # RMS error in coord positions, designed for SPI sim only but I guess it wouldn't be detrimental for crystal sim.
        self.random_waters = random_waters
        self.use_intensity_for_time=use_intensity_for_time



        self.ignore_deviations = False
        if self.positional_stdv == 0 and self.random_waters is None:
            self.ignore_deviations = True

        assert self.supercell_simulations <= self.num_supercells
        ## Parameters to be parsed by custom function because I cannot understand Bio.PDB.PDBParser's documentation.
        self.sym_rotations = []; self.sym_translations = [];   # Symmetry for each asymmetric unit simulated. If non-SPI, this is defining the supercell.
        self.cell_dim = None   # unit cell basis vector lengths. 
        self.cell_angles = None       
        self.parse_data_from_pdb() # All asymmetric units in unit cell

        if self.gromacs_config_file: 
            print("!!Do not use for scattering!!")
            self.cell_dim = [0,0,0]
            self.cell_angles = [90,90,90]

        if not include_symmetries:
            # Ignore parsed symmetries and use a single asymmetric unit per cell.
            self.sym_rotations = []; self.sym_translations = []; 
            self.add_symmetry_to_cells(np.identity(3),np.zeros((3)),"(X,Y,Z)")


        assert len(self.sym_rotations)!=0 and len(self.sym_translations)!=0
        
        if not self.ignore_deviations and self.use_bfactors:
            print("Warning: Using both deviations and B factors")

        self.supercell_dim = self.cell_dim*supercell_scale 


        ## Dictionary for going from pdb to ac4dc names.
        # different names
        PDB_to_AC4DC_dict = dict(
            #NA = "Sodion", CL = "Chloride",
            NA = "Na", CL = "Cl", CU="Cu",FE="Fe",CLA="CL",SOD="NA",CAL="Ca",CA="Ca",ZN="Zn"
        )
        # Same names
        for elem in ["H","He","C","N","O","P","S","Gd","I"]:  # pdb names # TODO automate this...
            PDB_to_AC4DC_dict[elem] = elem   # ac4dc name
        # Modify the values in the dictionary according to arguments.
        missing_species_names=[]
        missing_species_elements=[]
        for k,v in PDB_to_AC4DC_dict.items():
            # Light atom approximation. 
            if CNO_to_N:
                if v in ["C","O"]: 
                    v = "N"        
            if S_to_N:
                if v == "S":
                    v = "N"  
            if convert_excluded_elements_to_N and v not in allowed_atoms: 
                    v = "N"    
            if convert_excluded_elements_to_H and v not in allowed_atoms: 
                    v = "H"    

            
            PDB_to_AC4DC_dict[k] = v       
        if self.gromacs_config_file:
            parser = Custom_Gromacs_Parser()
        elif self.cif_file:
            parser = MMCIFParser()
            assert False, "cif formats not supported"
        else:
            # Get structure using Bio.PDB's parser
            parser=copy.deepcopy(xPDBParser)
        structure_id = os.path.basename(self.struct_file_path)
        try:
            with suppress_stdout_stderr():
                structure = parser.get_structure(structure_id, self.struct_file_path)  
        except Exception as e:
            print(f"error with {self.struct_file_path}")
            raise e
        self.num_input_file_atoms=len(list(structure.get_atoms()))
        # NOTE if modify here, need to modify in pdb_to_AC4DC_dict 
        for atom in structure.get_atoms(): # XXX Patch
            if atom.name in ("NA","CL","ZN"):
                atom.element = atom.name  
        # Get those cheeky charge clusters
        species_dict = {}
        pdb_atom_eles = []
        pdb_atoms_ignored = ""
        #TODO change to just reading data folder and only excluding atoms that are specified, rather than requiring user to pass in all allowed atoms.       
        for ac4dc_atom in allowed_atoms:
            if ac4dc_atom + "_fast" in allowed_atoms:
                #TODO read data folder to resolve ambiguity.
                raise Exception("Ambiguity, both "+ac4dc_atom+" and "+ac4dc_atom+"_fast were given in allowed_atoms.") 
            if ac4dc_atom + "_faster" in allowed_atoms:
                #TODO read data folder to resolve ambiguity.
                raise Exception("Ambiguity, both "+ac4dc_atom+" and "+ac4dc_atom+"_faster were given in allowed_atoms.")            
        charge_idxes={}

        # with open(self.struct_file_path) as f:
        #     for line in f:
        #         print(line)
        # print("STRUCTURE FILE PATH:", self.struct_file_path)
        # print("ATOM NAMES:", [a.get_name() for a in structure.get_atoms()])
        for i, atom in enumerate(structure.get_atoms()):
            if ignore_water_H and atom.element=="H" and atom.get_parent() is not None and atom.get_parent().get_resname()=="HOH":
                continue        
    
            special_convert_dict={"FE":"FE","CLA":"CL","SOD":"NA","ZN":"ZN","CAL":"CA"} # changes here change in gromacs parser
            if atom.get_name() in special_convert_dict:
                atom.element = special_convert_dict[atom.get_name()]            
            if (atom.get_name() == "NA" and atom.get_parent() is not None 
            and atom.get_parent().get_resname()=="HEM"):
                atom.element="N"
            # Get ze data
            R = atom.get_vector()
            ac4dc_name = PDB_to_AC4DC_dict.get(atom.element)
            # Pop into our desired format
            if ac4dc_name == None:
                pdb_atoms_ignored += atom.element + " "
                continue
            # TODO make this not suck.
            original_ac4dc_name=ac4dc_name
            if ac4dc_name not in allowed_atoms: 
                ac4dc_name+= "_fast"
                if ac4dc_name not in allowed_atoms: 
                    ac4dc_name +="er"
                    if ac4dc_name not in allowed_atoms:
                        ac4dc_name = original_ac4dc_name
                        ####
                        if not allow_skip_species: 
                            assert False, f"{ac4dc_name} was not in expected elements: {allowed_atoms}" 
                        else:
                            missing_species_names.append(ac4dc_name)
                            missing_species_elements.append(atom.element)
                            allowed_atoms.append(ac4dc_name)
                        ####
            if ac4dc_name not in species_dict.keys():
                PDB_to_AC4DC_dict[atom.element]=ac4dc_name
                species_dict[ac4dc_name] = Atomic_Species(ac4dc_name,self) 
                charge_idxes[ac4dc_name]=[]
                pdb_atom_eles.append(atom.element)
            species_dict[ac4dc_name].add_atom(atom.get_serial_number(),R,atom.get_bfactor())
            charge_idxes[ac4dc_name].append(i)
        if charge_states is not None:
            for ac4dc_name in species_dict:
                species_dict[ac4dc_name].set_charge_states(charge_states[charge_idxes[ac4dc_name]])
            

        ac4dc_atoms_ignored = ""
        for string in allowed_atoms:
            if string not in species_dict.keys():
                ac4dc_atoms_ignored += string + " "
                continue
        print("The following atoms will be considered:\nAC4DC names ; pdb names (num / asu)")
        for ac,pd in zip([PDB_to_AC4DC_dict[x] for x in pdb_atom_eles if x not in missing_species_elements], pdb_atom_eles): 
            print(f"{ac:<10} {pd:>10} ({species_dict[ac].get_num_atoms()})")
        if len(missing_species_names)>0:
            print("WARNING: The following atoms are present in the structure but ignored:\n(AC4DC names ; pdb names)")
            for p,a in zip([PDB_to_AC4DC_dict[x] for x in missing_species_elements], missing_species_elements): 
                print(f"{p:<10} {a:>10}")
        if pdb_atoms_ignored != "":
            print("The following pdb atoms were found but ignored:",list(set(pdb_atoms_ignored)))
        if ac4dc_atoms_ignored != "":
            print("The following atoms were allowed but not found:",ac4dc_atoms_ignored)
        
        self.missing_species_dict = {k:v for k,v in species_dict.items() if k in missing_species_names}         
        self.species_dict = {k:v for k,v in species_dict.items() if k not in missing_species_names} 

        print("Number of atoms no symm",self.num_atoms_no_symm())

        #self.set_stochastic_positions(first_call=True)
    def disable_pos_deviations(self):
        self.positional_stdv=0
        self.ignore_deviations=True

    def num_atoms_no_symm(self):
        return np.sum([len(self.species_dict[k].coords) for k in self.species_dict])
    def set_stochastic_positions(self,q):
        first_call=True
        if self.stochastic_positions_set:
            self.stochastic_positions_set=True
            first_call = False
        if self.random_waters is not None:
            if first_call:
                self.random_water_start_idx = self.num_atoms_no_symm()
                self.num_non_water_oxygens = 0 if "O" not in self.species_dict else len(self.species_dict["O"].coords)
            self.reinitialize_random_waters()
            if first_call:
                print(f"Added {self.random_waters} O atoms to random coordinates. Structure now has {self.num_atoms_no_symm()} atoms (asu).")
        else:
            print(f"Structure has {self.num_atoms_no_symm()} atoms (asu).")
        if first_call:
            print(f"Adding stdv of {self.positional_stdv*ang_per_bohr} Angstroms.")
        for species in self.species_dict.values():
            species.set_coord_deviation(q)
        print("coord deviations set")

    def reinitialize_random_waters(self):
        if "O" not in self.species_dict:
            return
        self.species_dict["O"].coords = self.species_dict["O"].coords[:self.num_non_water_oxygens]
        if DEBUG or DEBUG_WATER:
            print(f"Placing {self.random_waters} O atoms in random positions")
        for _ in range(self.random_waters):
            random_water_bfactor=100
            self.species_dict["O"].add_atom("WATER",Bio_Vect((np.random.rand(3)-0.5)*1e3),random_water_bfactor) # make water b factor high instead?
    def set_ff_calculator(self,ff_calculator : Plotter):
        self.ff_calculator = ff_calculator                  
    
    def add_symmetry_to_cells(self,symmetry_factor,symmetry_translation,symmetry_label="Symmetry Operation:"):
        '''
        Takes in the symmetry factor and translation of the unit cell, and translates them across each unit cell
        '''
        symmetry_translation/= ang_per_bohr

        if (self.cell_packing == "SC" or self.cell_packing == "triclinic")\
            and (not self.cell_angles[0] == self.cell_angles[1] == self.cell_angles[2]==90):
            print("Warning: Overriding with triclinic basis")
            self.cell_packing="triclinic"



        # if self.cell_packing == "triclinic":
        #     symmetry_factor = symmetry_factor @ get_triclinic_basis(self.cell_angles) #get_triclinic_basis(self.cell_angles)@symmetry_factor
        # Simple cubic packing. (actually it's all rectangular prisms, TODO)
        if self.cell_packing == "SC" or self.cell_packing == "triclinic":  
            if self.cell_packing == "SC":
                assert self.cell_angles[0] == self.cell_angles[1] == self.cell_angles[2]==90

            self.num_cells = self.supercell_scale**3
            # Generate the coordinates of the cube (performance: defining scaling matrix/cube coords outside of here would be more efficient). But only runs once  \_( '_')_/ ¯\_(ツ)_/¯.
            x, y, z= np.meshgrid(np.arange(0, self.supercell_scale), np.arange(0, self.supercell_scale), np.arange(0, self.supercell_scale))
            cube_coords = np.stack([ x.flatten(), y.flatten(), z.flatten()], axis = -1)
            # Construct the cube by adding the "unit cell translation" to the symmetry translation. 
            # (e.g. if crystal is centred on origin, in fractional crystallographic coordinates the index of the unit cell is equal to the unit cell's translation from the origin.) 
            if DEBUG or DEBUG_MODERATE:
                print_thing = np.array(["][x]   [","][y] + [","][z]   ["],dtype=object)
                print(symmetry_label,"\n",np.c_[symmetry_factor, print_thing ,symmetry_translation],sep="")
          
            if self.cell_packing == "triclinic":
                a = get_triclinic_basis(self.cell_angles)
            for coord in cube_coords:
                if self.cell_packing == "triclinic":
                    coord_scaled = (coord*self.cell_dim)@a
                    assert(coord.shape == (3,))
                    translation = coord_scaled + symmetry_translation

                #print(coord)
                else:
                    translation = coord*self.cell_dim + symmetry_translation
                self.sym_translations.append(translation)
                self.sym_rotations.append(symmetry_factor)
                 
        else:
            raise Exception("Lacking implementation")
            
    def parse_data_from_pdb(target):
        '''
        Parses the symmetry transformations from the pdb file into ndarrays and adds each to the crystal.  
        No unit operations are performed.
        Also returns the parsed matrices. 
        '''
        at_SYMOP = False
        at_symmetry_xformations = False
        sym_factor = np.zeros((3,3))
        sym_trans = np.zeros((3))
        sym_labels = []
        sym_mtces_parsed = []        
        with open(target.struct_file_path) as pdb_file:
            for line in pdb_file:
                line = line.strip()
                # End data - Check if left section of data; flagged by the line containing solely "REMARK ###". 
                if len(line) <= 10: 
                    if at_SYMOP:
                        at_SYMOP = False
                    if at_symmetry_xformations: 
                        at_symmetry_xformations = False
                if line[0:4] == "ATOM":
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
                        sym_mtces_parsed.append([sym_factor.copy(),sym_trans.copy()])
                        

                ## inline data
                if line[0:6] == "CRYST1":
                    print(target.cell_dim)
                    if target.cell_dim is not None:
                        raise Exception("Cannot handle multiple crystal entries")
                    entries = line.split()[1:]
                    target.cell_dim = [float(a) for a in entries[0:3]]
                    target.cell_dim = np.array(target.cell_dim)/ang_per_bohr 
                    target.cell_angles = [float(a) for a in entries[3:6]]
                    target.symmetry = entries[6:-1]


                ## Start data - Check if this line marks the next as the beginning of a desired data section.                    
                # Symmetry operators
                if "NNNMMM   OPERATOR" in line:
                    at_SYMOP = True
                    if DEBUG or DEBUG_MODERATE:
                        print("Parsing symmetry operations...")
                # symmetry matrices
                if line == "REMARK 290 RELATED MOLECULES.":
                    at_symmetry_xformations = True
        for sym_factor, sym_trans in sym_mtces_parsed:
            target.add_symmetry_to_cells(sym_factor.copy(),sym_trans.copy(),sym_labels.pop(0))
        return sym_mtces_parsed
    
    unique_points=None
    def reset_unique_points(self):
        self.unique_points = None
    def get_sym_xfmed_point_unique_mask(self,R,symmetry_index):
        i = symmetry_index
        points = self.get_sym_xfmed_point(R,i)
        if self.unique_points is None:
            self.unique_points = np.array([],dtype=points.dtype).reshape(0,3)
   
        is_unique_check = (self.unique_points[:, None] == points).all(axis=2).any(axis=0) ==False
        points = points[is_unique_check]
        # for i, point in enumerate(points):
        #     if point in self.unique_points:
        #         del points[i]
        self.unique_points = np.append(self.unique_points,points,axis=0)
                
        assert np.unique(self.unique_points,axis=0).shape==self.unique_points.shape, f"{np.unique(self.unique_points,axis=0).shape} {self.unique_points.shape}"
        return is_unique_check
        

    def get_sym_xfmed_point(self,R,symmetry_index):
        '''
        R = 1D or 2D array of shape (3,N)
        '''
        i = symmetry_index
        # Transpose because we've stored the vectors in the 0th axis. 
        return (self.sym_rotations[i] @ R.T).T+ self.sym_translations[i]   # dim = [xyz,xyz] x [xyz,N] or dim = [xyz,xyz] x [xyz,1]
            
    def save_structure_by_reference(self,dir="targets",tag="constructed_struct",allow_partial_occupancy=False):
        '''
        Version of save_structure with less misunderstandings
        '''
        #Load it
        parser=copy.deepcopy(xPDBParser)
        reference_structure_id = os.path.basename(self.struct_file_path)
        reference_structure = parser.get_structure(reference_structure_id, self.struct_file_path)    

        #Initialise structure to build
        structure = xStructureBuilder()
        structure.init_structure(tag)
        structure.init_model("M")
        structure.init_seg("")

        # Add atoms
        serial_number = 1
        chainIDs = ascii_uppercase+ascii_lowercase+digits
        num_chains = len(self.sym_rotations)
        one_chain_per_unit=False

        altloc_mode=False
        num_reference_chains=len(list(reference_structure.get_chains()))
        if not altloc_mode:
            if num_chains > len(chainIDs) or num_reference_chains>1:
                one_chain_per_unit = True
                num_chains = self.supercell_scale**2*num_reference_chains
                assert num_chains < len(chainIDs)
        asym_per_cell = len(self.sym_rotations)/self.supercell_scale**2
        residues_per_asym=len(list(reference_structure.get_residues()))
        if num_chains < len(ascii_uppercase):
            chainIDs = ascii_uppercase
        abs_chain_idx=-1
        chains_list = list(reference_structure.get_chains())
        for i in range(len(self.sym_rotations)):
            if  (i%asym_per_cell == 0):
                chain_last_resnum_dict={}

            unit_cell_idx=int(np.floor(i/asym_per_cell))
            if altloc_mode:
                altloc = ascii_letters[unit_cell_idx]
            else:
                altloc = ' '
            print(altloc)
            for c,chain in enumerate(chains_list):
                abs_chain_idx+=1
                #abs_chain_idx=num_reference_chains*asym_per_cell*c+i
                if one_chain_per_unit:
                    chain_idx=unit_cell_idx
                else:
                    # Share chain ids between unit cells
                    chain_idx=int(abs_chain_idx-unit_cell_idx*asym_per_cell*num_reference_chains)
                chainID=chainIDs[chain_idx]
                structure.init_chain(chainID) 
                
                for j, reference_residue in enumerate(chain.get_residues()):
                    hetflag, resseq, icode = reference_residue.get_id()
                    if chainID not in chain_last_resnum_dict:
                        chain_last_resnum_dict[chainID]=0
                    chain_last_resnum_dict[chainID]+=1
                    resseq=chain_last_resnum_dict[chainID]
                    assert resseq <1e4, f"Residue number {resseq} too high!"
                    print(resseq)

                    r_args = (hetflag,resseq,icode)
                    #get_resname()
                    #get_segid()
                    
                    residue_id = r_args
                    if residue_id not in structure.chain:
                        structure.init_residue(reference_residue.get_resname(),*r_args)
                    residue = structure.chain[residue_id] 
                    for reference_atom in reference_residue.get_atoms():
                        R = reference_atom.get_vector().get_array()/ang_per_bohr

                        coord = self.get_sym_xfmed_point(R,i)
                        
                        coord=tuple([c*ang_per_bohr for c in coord])

                        if not allow_partial_occupancy:
                            assert reference_atom.get_occupancy()==1, (reference_atom.get_occupancy(),reference_atom.full_id)
                        atom = PDB_Atom(name=reference_atom.get_name(), coord=coord, bfactor=reference_atom.get_bfactor(), occupancy=reference_atom.get_occupancy(), 
                                                altloc=altloc, fullname=reference_atom.get_fullname(), serial_number=serial_number,element=reference_atom.element)
                        disordered_atom=None
                        for a in residue.get_atoms():
                            if a is None:
                                continue
                            if reference_atom.get_name()==a.get_name():
                                disordered_atom=a
                        if disordered_atom is None:
                            disordered_atom=DisorderedAtom(atom.get_name())
                            residue.add(disordered_atom)  
                        disordered_atom.disordered_add(atom)
                        serial_number+=1

        
        # Save it
        io=xPDBIO()
        io.set_structure(structure.get_structure())  # StructureBuilder object is not Structure object
        fname = path.basename(self.struct_file_path)[:-4]+f"_{tag}.pdb"
        io.save(dir+'/'+fname)


        a,b,c = self.cell_dim*ang_per_bohr 
        header=f"""HEADER    full_struct                          date   name              
TITLE     FULL_STRUCT    

REMARK 290 SYMMETRY OPERATORS FOR SPACE GROUP: P 1                        
REMARK 290                                                                      
REMARK 290      SYMOP   SYMMETRY                                                
REMARK 290     NNNMMM   OPERATOR                                                
REMARK 290       1555   X,Y,Z                                                                                            
REMARK 290                                                                      
REMARK 290     WHERE NNN -> OPERATOR NUMBER                                     
REMARK 290           MMM -> TRANSLATION VECTOR                                  
REMARK 290                                                                      
REMARK 290 CRYSTALLOGRAPHIC SYMMETRY TRANSFORMATIONS                            
REMARK 290 THE FOLLOWING TRANSFORMATIONS OPERATE ON THE ATOM/HETATM             
REMARK 290 RECORDS IN THIS ENTRY TO PRODUCE CRYSTALLOGRAPHICALLY                
REMARK 290 RELATED MOLECULES.                                                   
REMARK 290   SMTRY1   1  1.000000  0.000000  0.000000        0.00000            
REMARK 290   SMTRY2   1  0.000000  1.000000  0.000000        0.00000            
REMARK 290   SMTRY3   1  0.000000  0.000000  1.000000        0.00000                        
REMARK 290   
CRYST1   {a*self.supercell_scale:.3f}   {b*self.supercell_scale:.3f}   {c*self.supercell_scale:.3f}  {self.cell_angles[0]:.2f} {self.cell_angles[1]:.2f}  {self.cell_angles[2]:.2f}   P 1    1"""

        with open(dir+'/'+fname, 'r+') as f:
            content = f.read()
            
            print(len(header.split('\n')))
            f.seek(0, 0)
            f.write(header.rstrip('\r\n') + '\n' + content)

        print(f"Saved structure to {dir+'/'+fname}")

    
    def save_structure(self,dir="targets",custom_residue_name=None,tag="constructed_struct",chain_name="C"):
        '''
        Saves the full structure in a pdb file format for use with Solvate1.0   
        Stdv not included.
        I have not tested if using it for over 1e5 atoms (when xpdb alters things so it doesn't break) works when plugged into SOLVATE
        TODO need to add boilerplate (replace HETATM with 'ATOM  ', add symmetries and cell dimensions to start of file.)
        '''
        #Load it
        # parser=PDBParser(PERMISSIVE=1)
        # structure_id = os.path.basename(self.struct_file_path)
        # structure = parser.get_structure(structure_id, self.struct_file_path)       

        #Initialsie structure
        structure = xStructureBuilder()
        structure.init_structure(tag)
        structure.init_model("M")
        structure.init_seg("")

        # Add atoms
        serial_number = 1
        for i in range(len(self.sym_rotations)):
            residue_name = "L"+str(i); r_args = (" ",i+1,"r")
            if custom_residue_name is not None:
                residue_name = custom_residue_name
            print(structure.chain.child_dict)
            #residue_id = (r_args[0]+"_"+residue_name,r_args[1],r_args[2])  # use if set r_args[0] to "H"
            residue_id = r_args
            structure.init_residue(residue_name,*r_args)
            print(structure.chain.child_dict)
            residue = structure.chain[residue_id] 
            for species in self.species_dict.values():
                points = np.array(species.coords)                    
                coord_list = self.get_sym_xfmed_point(points,i).tolist()
                for j, coord in enumerate(coord_list):
                    coord=tuple([c*ang_per_bohr for c in coord])
                    name = ' '+species.name+str(j)+' '
                    residue.add(PDB_Atom(name= name, coord=coord, bfactor=0., occupancy=1., altloc=' ', fullname=name, serial_number=serial_number,element=species.name))                
                    serial_number+=1
        # TODO boiler plate
        # structure.set_symmetry("P 1", "1")
        
        # Save it
        io=xPDBIO()
        io.set_structure(structure.get_structure())  # StructureBuilder object is not Structure object
        fname = path.basename(self.struct_file_path)[:-4]+f"_{tag}.pdb"
        io.save(dir+"/"+fname)     

    def plot_me(self,max_points = 100000,water_index = None,**layout_kwargs):
        if water_index != None:
            assert self.supercell_scale == 1 and len(self.sym_rotations) == 1, "Plotting water with an added intra-cell or crystal symmetry is not supported."
        '''
        plots the symmetry points. Origin is the centre of the base asymmetric unit (not the centre of a unit cell).
        RMS error in positions ('positional_stdv') is not accounted for
        ''' 
        plt.close()
        view_width = 1000
        view_height = 800

        if water_index is None and self.random_waters is not None:
            water_index = self.random_water_start_idx # TODO make work with symmetries
        num_atoms_avail = 0
        for species in self.species_dict.values():
             num_atoms_avail += len(species.coords)
        assert water_index is None or water_index < num_atoms_avail
        num_test_points = max(1,int(max_points/(len(self.sym_translations))))
        num_test_points = min(num_atoms_avail,num_test_points)         
        test_points = np.empty((num_test_points,3))
        qualifier = "sample of"
        if num_atoms_avail == len(test_points):
            qualifier = "all" 
        print("Number of atoms in supercell:",num_atoms_avail*len(self.sym_translations))
        print("Generating plot with",qualifier,len(test_points),"atoms per asymmetric unit (plotting "+str(len(test_points)*len(self.sym_translations))+" in total)")
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

        plot_coords = np.array(plot_coords)*ang_per_bohr # convert to angstrom
        raise_non_unique_exception = False
        if np.unique(plot_coords,axis=0).shape != plot_coords.shape:
            a, unique_indices = np.unique(plot_coords,axis=0,return_index = True)
            num_non_unique = len(plot_coords)-len(unique_indices)
            print("WARNING", num_non_unique, "non-unique coords found:")
            if not DELETENONUNIQUE:
                raise_non_unique_exception = True
            #plot_coords = np.delete(plot_coords, unique_indices,axis=1)  # delete anyway.            
            #plot_coords = plot_coords[unique_indices]

        #Colors - which to highlight (root atom is first atom of cell)
        c_first_root_atom = False # red
        c_first_unit = True # pink
        c_first_cell = True # aquamarine
        c_all_root_atoms = False # yellow

        top_atom_index = 0 # index of atom that has highest z
        max_height = -np.inf
        for i,coord in enumerate(test_points):
            if coord[2] > max_height:
                max_height = coord[2]
                top_atom_index = i 
        if water_index == None:
            alpha = None
            color = np.empty(len(plot_coords),dtype=object) 
            color.fill('green')
            if c_first_cell:
                for i in range(int(len(self.sym_rotations)/self.num_cells)):           
                    color[i*self.num_cells*len(test_points):i*self.num_cells*len(test_points)+len(test_points)] = 'aquamarine' #'c'  # atoms in one unit cell  
            if c_first_unit:
                color[:len(test_points)] = 'pink'  # atoms in one asym unit       
            #Cursed but it works. 
           
            for i in range(self.num_cells):
                if c_all_root_atoms:
                    for j in range(int(len(self.sym_rotations)/self.num_cells)):
                        color[i*int(len(test_points)*len(self.sym_rotations)/self.num_cells) + j*len(test_points) + top_atom_index] = 'yellow'      # same atom in every asym unit.  
                if c_first_root_atom:
                    color[i*len(test_points)+top_atom_index] = 'red'      # same atom in every same unit cell         
        else:
            kernel_color = np.empty(water_index+1,dtype=object);# kernel_alpha =  np.empty(water_index+1,dtype=np.double)
            bg_color = np.empty(len(plot_coords)-(water_index+1),dtype=object); #bg_alpha = np.empty(len(plot_coords)-(water_index+1),dtype=np.double)
            kernel_color.fill('rgba(92,169,4,1)'); #kernel_alpha.fill(1) 
            bg_color.fill('rgba(0,30,255,0.2)');  # water atom colours
            color = np.concatenate((kernel_color,bg_color))
            #alpha = np.concatenate((kernel_alpha,bg_alpha))
        

            #color[color!='y'] = '#0f0f0f00'   # example to see just one atom type
        x_range = np.array([np.min(plot_coords[:,0]),np.max(plot_coords[:,0])])
        y_range = np.array([np.min(plot_coords[:,1]),np.max(plot_coords[:,1])])
        z_range = np.array([np.min(plot_coords[:,2]),np.max(plot_coords[:,2])])
        print("---")
        print("x_min/max:",x_range[0],x_range[1])
        print("y_min/max:",y_range[0],y_range[1])
        print("z_min/max:",z_range[0],z_range[1])
        print("---")
        max_len = 0
        min_len = np.inf
        for elem in x_range,y_range,z_range:
            max_len = max(max_len,abs(elem[1]-elem[0]))      
            min_len = min(min_len,abs(elem[1]-elem[0]))      

        titles = [{"text": ax+'$ [\AA]$'} for ax in ['x','y','z']]
        xaxis,yaxis,zaxis = [{'title': title} for title in titles]
        ##
        custom_axes = False #(i.e. we won't use the axis variables above)
        ranges = np.array([np.max(plot_coords[:,i]) - np.min(plot_coords[:,i]) for i in range(3)])
        if custom_axes:
            # data for custom axis centred in the space
            avg_range = np.sum(ranges)/3
            step_size = max(1, 2.5*10**(np.log10(avg_range/10) - ((np.log10(avg_range/10)-1)%1)))
            # Translate to be centred on axis ticks
            plot_coords -= np.array([np.mean(plot_coords[:,i]) for i in range(3)])        
            centre_x = centre_y = centre_z = 0
            #centre_x = np.mean(plot_coords[:,0]); centre_y = np.mean(plot_coords[:,1]); centre_z = np.mean(plot_coords[:,2])

            #centre_x -= centre_x%step_size; centre_y -= centre_y%step_size; centre_z -= centre_z%step_size 
        origin_on_corner = True  
        # asym unit
        # size = 5.5*max(1,500/(min_len+max_len))   
        # angular_aperture = np.pi*0.6 
        # dot_lw = 1
        # unit cell
        # size = 7.5*max(1,500/(min_len+max_len))   
        # angular_aperture = np.pi*0.6
        # dot_lw = 0.1 # unit
        #2x2x2 Crystal
        # size = 10*max(1,500/(min_len+max_len))   #solvated crystal
        # angular_aperture = np.pi*0.7 # Solvated crystal        
        # dot_lw = 0 # cryst    

        # tetrapeptide:
        # size = 3*max(1,500/(min_len+max_len))
        # angular_aperture = np.pi*0.6
        # dot_lw = 0.1 

        # non-camera perspective (dots all same):
        size = 1.5*max(1,500/(min_len+max_len))
        dot_lw = 0
        if origin_on_corner:
            max_num_ticks = 7
            minima = np.array([np.min(plot_coords[:,i]) for i in range(3)])
            
            #Translate origin to corner.
            plot_coords -= minima  
            minima = np.array([np.min(plot_coords[:,i]) for i in range(3)])
            maxima = np.array([np.max(plot_coords[:,i]) for i in range(3)])
            tick_spacing = 10
            addendums = [dict(tickfont=dict(size=15), showbackground=False,mirror="all", nticks=1+max(1,round(m*1.4/tick_spacing)), range=[0,m*1.4]) for m in maxima] #nticks=max(2,round(m/max(maxima)*max_num_ticks))
            addendums[1]["range"] = [0,maxima[1]] #pull back wall in
            addendums[1]["nticks"] = 1+max(1,round(maxima[1]/tick_spacing))
            for i, axis in enumerate([xaxis,yaxis,zaxis]):
                axis |= addendums[i]
            
            # INCREDIBLY LAZY implementation to provide some sense of distance in png snapshots.
            hack_snapshots = False
            if hack_snapshots:
                # Calculate size relative to view from camera that is capturing full width of structure.
                # Hacky as I am assuming camera is at a certain position.
                r = size  # some arbitrary 'radius' of visible atom that is equivalent to the largest possible pixel width. 
                #width = view_width/r
                # a guess (increase to make fall off more slowly)

                target_length = np.sqrt(np.sum([range**2 for range in ranges]))
                distance_to_full_capture = np.tan(angular_aperture/2)*(target_length/2)
                #psi = -np.pi/4; thet=np.pi/4 #Camera position relative to target's centre
                psi = -np.pi/2; thet=np.pi/6 # cryst alternate angle
                camera_pos = [distance_to_full_capture*np.sin(thet)*np.cos(psi),distance_to_full_capture*np.sin(thet)*np.sin(psi),distance_to_full_capture*np.cos(thet)]
                max_atom_angle = np.arctan(r/(target_length/2))
                atom_distances = [np.sqrt(np.sum([(val-camera_pos[i])**2 for i, val in enumerate(coord)])) for coord in plot_coords]
                #r*angle/max_angle
                size = [r*np.arctan(r/(distance_to_full_capture + d))/max_atom_angle for d in atom_distances]
        #size = max(1,500/(min_len+max_len))  


        scatter_points = go.Scatter3d(
            x=plot_coords[:,0], 
            y=plot_coords[:,1], 
            z=plot_coords[:,2], 
            marker=go.scatter3d.Marker(color=color,size=size,line=dict(width=dot_lw,color='black'),opacity=1), 
            mode='markers',
        )
        fig=go.Figure(data=scatter_points)
        # def sphere(x, y, z, radius, resolution=4):
        #     """Return the coordinates for plotting a sphere centered at (x,y,z)"""
        #     u, v = np.mgrid[0:2*np.pi:resolution*2j, 0:np.pi:resolution*1j]
        #     X = radius * np.cos(u)*np.sin(v) + x
        #     Y = radius * np.sin(u)*np.sin(v) + y
        #     Z = radius * np.cos(v) + z
        #     return (X, Y, Z)
        # scatter_spheres = []
        # for point in plot_coords:
        #     (x_pns_surface, y_pns_surface, z_pns_surface) = sphere(*point,max(1,500/(min_len+max_len)))
        #     scatter_spheres.append(go.Surface(x=x_pns_surface, y=y_pns_surface, z=z_pns_surface, opacity=0.5))
        #fig=go.Figure(data=scatter_spheres)
        # Setup 3D scene stuff
        aspect = np.empty(3)
        for i, elem in enumerate((x_range*1.4,y_range,z_range*1.4)):
            aspect[i] = abs(elem[1]-elem[0])/max_len    
        zoom = 1.8 # inital zoom      
        fig.update_layout(
            margin=dict(l=20,r=20,t=20,b=20),
            width=view_width, height=view_height,
            scene = dict(
                xaxis = xaxis,
                yaxis = yaxis,
                zaxis = zaxis,
                aspectratio = dict(x=aspect[0]*zoom,y=aspect[1]*zoom,z=aspect[2]*zoom),
                #camera=dict(eye=dict(x=1,y=0,z=0.6))
            ),
        ) 
        if custom_axes:
            # Remove default axis
            axis_args = dict(show_grid = False, zeroline = False, showticklabels = False, title = dict(text = ""))
            fig.update_layout(scene = dict(zaxis=axis_args,yaxis=axis_args,xaxis=axis_args))      
            # Add custom axis centred in the space
            x_tickvals = np.append(
                np.flip(np.arange(centre_x,np.min(plot_coords[:,0])+abs(np.min(plot_coords[:,0]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_x,np.max(plot_coords[:,0])-abs(np.max(plot_coords[:,0]))%step_size+step_size*3/4,step_size)
            )
            y_tickvals = np.append(
                np.flip(np.arange(centre_y,np.min(plot_coords[:,1])+abs(np.min(plot_coords[:,1]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_y,np.max(plot_coords[:,1])-abs(np.max(plot_coords[:,1]))%step_size+step_size*3/4,step_size)
            )
            z_tickvals = np.append(
                np.flip(np.arange(centre_z,np.min(plot_coords[:,2])+abs(np.min(plot_coords[:,2]))%step_size-step_size*3/4,-step_size)),
                np.arange(centre_z,np.max(plot_coords[:,2])-abs(np.max(plot_coords[:,2]))%step_size+step_size*3/4,step_size)
            )
            x_tickvals = np.round(x_tickvals,0)
            y_tickvals = np.round(y_tickvals,0)
            z_tickvals = np.round(z_tickvals,0)
            print(np.flip(np.arange(centre_y,np.min(plot_coords[:,1])+abs(np.min(plot_coords[:,1]))%step_size,-step_size)))
            line_width = 10
            marker_size = 3
            fontsize = 20
            xaxis_line =go.Scatter3d(
                            x = x_tickvals,
                            y = (centre_y,)*len(x_tickvals),
                            z = (centre_z,)*len(x_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=x_tickvals,
                            textfont=dict(size=fontsize),
                            )
            yaxis_line =go.Scatter3d(
                            y = y_tickvals,
                            z = (centre_z,)*len(y_tickvals),
                            x = (centre_x,)*len(y_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=y_tickvals,
                            textfont=dict(size=fontsize),
                            )
            zaxis_line =go.Scatter3d(
                            z = z_tickvals,
                            x = (centre_x,)*len(z_tickvals),
                            y = (centre_y,)*len(z_tickvals),
                            mode = "lines+markers+text",
                            marker = dict(size=marker_size),
                            line = dict(width = line_width),
                            text=z_tickvals,
                            textfont=dict(size=fontsize),
                            )                
            centre_ball = go.Scatter3d(x = (centre_x,centre_x), 
                                    y = (centre_y,centre_y), 
                                    z = (centre_z,centre_z), 
                                    mode = "markers", 
                                    hoverinfo = "skip", 
                                    marker = dict(size = 8))
            fig.add_traces([xaxis_line,yaxis_line,zaxis_line,centre_ball])

        # Update the kwargs separately so that we can call arguments used above without conflict.
        fig.update_layout(**layout_kwargs)
        fig.show()       
        #fig.write_html("Structure.html")

        if raise_non_unique_exception:
            raise Exception(str(num_non_unique) + " non-unique coords found:")


class Atomic_Species():
    def __init__(self,name,crystal: Crystal):
        self.name = name 
        self.crystal = crystal 
        self.coords = []  # coord of each atom in species in asymmetric unit
        self.serial_numbers = [] # Corresponding serial number of each atom 
        self.B_factors = [] 
        self.debye_waller_factor=None

        self.charge_states=None # Charge states for all times and atom indices
        self.ff_by_occupancy_and_time=None
        # if self.charge_states is not None:
        #     assert len(charge_states) == len() 

    def set_charge_states(self,charge_states):
        self.charge_states=charge_states
    def Z(self):
        return ATOMNO[self.name]

    def add_atom(self,serial_number,vector,B_factor):
        '''
        This function adds an atom to the asymmetric unit of the crystal. 
        We do not store additional coordinates, instead storing the symmetries, and an array of atomic states corresponding to each atom, for each symmetry. (so num symmetries * num atoms added)
        '''           
        self.serial_numbers.append(serial_number)
        self.coords.append(vector.get_array()/ang_per_bohr)
        if not self.crystal.zero_bfactors:
            self.B_factors.append(B_factor/(ang_per_bohr**2))
    def set_debye_waller(self,q):
        if self.crystal.zero_bfactors:
            return 1
        assert len(self.B_factors) == len(self.coords)
        factor = np.array(self.B_factors)/(16*np.pi**2) 
        if len(q.shape)==1:
            self.debye_waller_factor = np.exp(-np.square(q)[None,...]*factor[:,None])
        elif len(q.shape)==2: # better way to do this...?
            self.debye_waller_factor = np.exp(-np.square(q)[None,...]*factor[:,None,None])

        #print(np.array(self.B_factors))
        #print(self.debye_waller.shape)
        
    def set_stochastic_electronic_states(self,q_arr = None):
        '''
        We set a state for each atom, including symmetries, so that different q applied to the same atom at the same time corresponds to the same state. 
        When we are finding the atomic form factors, we call get_atomic_form_factors() on these states.
        The dimensionless nature of the model forces us to make the dubious approximation that an atom's state is independent of its prior states.
        With a hybrid molecular dynamics model informed by AC4DC, the nuclei's states could potentially be tracked properly through time, and this function would be replaced
        by a call to the data of the atomic nuclei's states.
        '''
        if DEBUG or DEBUG_MODERATE:
            print("Creating time-varying states for atom "+self.name+" from plasma simulation's data")
        self.times_used = self.crystal.ff_calculator.get_times_SCATTER()
        #print(f"Snapshot times:",self.times_used)
        num_atoms = self.get_num_atoms()
        
        if ( self.B_factors and (len(self.B_factors)!=len(self.coords))) \
        or ((not self.B_factors and not self.crystal.zero_bfactors) and  (self.num_atoms_on_coord_deviation != num_atoms)):
            raise Exception("num atoms was not same on set_stochastic_electronic_states call as when set by set_coord_deviation")
        
        
        if self.crystal.is_damaged:
            # Initialise an array that tracks the form factors of individual atoms.
            orb_occs_shape = (num_atoms,len(self.times_used))  # [num atoms,times]
            # TODO instead of storing lists, replace with indices and a list with corresponding states. Also use index to get ff rather than orbocc list
            self.orb_occs= np.empty(orb_occs_shape,dtype=list)    # self.orb_occs[i] is an array of states corresponding to each time. We make the necessary approximation that an atom's state is independent of its prior states. This approximation is dubious at low unit cell numbers, but at higher numbers, because the contribution from an atom at the same relative cell coordinate and state as another atom will be equivalent, we get the same outcome so long as the probability distribution of states is representative of the actual distribution of states. i.e. tracking state history at the same global position is redundant at high unit cell counts where we can expect the distribution of states at a given coordinate to have a low deviation between species.  
            
            if self.charge_states is None:
                for idx in range(num_atoms): # XXX Super slow
                    seed = None
                    if SEEDED:
                        seed = idx
                    self.orb_occs[idx],_,self.orb_occ_dict = self.crystal.ff_calculator.random_state_snapshots(self.name,seed) 
            else: 
                # NOTE: Not stochastic
                if self.ff_by_occupancy_and_time is not None:
                    # Already set. 
                    return 
                assert q_arr is not None
                self.atom_occupancies=self.Z()-self.charge_states
                total_occupancies:dict[int,list]={}
                config_strings= [occ_str for occ_str in self.crystal.ff_calculator.statedict[self.name]]
                orb_occs = [parse_elecs_from_latex(occ_str) for occ_str in config_strings]
                for config_str,occ in zip(config_strings,orb_occs):
                    total_occupancy = np.sum(list(occ.values()))
                    if total_occupancy not in total_occupancies:
                        total_occupancies[total_occupancy]=[]
                    total_occupancies[total_occupancy].append(config_str)
                self.occupancy_indices={total_occ:i for i, total_occ in 
                                    enumerate(sorted(list(total_occupancies.keys())))}
                self.ff_by_occupancy_and_time = np.zeros(
                    shape=(
                        len(total_occupancies), 
                        len(self.times_used),
                        *q_arr.shape))
                t_idx = np.searchsorted(self.crystal.ff_calculator.timeData,self.times_used)     
                


                #total_density=np.sum(self.crystal.ff_calculator.boundData[self.name][0,:])
                occ_dict=self.crystal.ff_calculator.get_occ_dict(self.name)
                got_ff_once=False
                assert len(total_occupancies)>0
                self.unused_occupancy_indices=[]
                for occupancy,configs in total_occupancies.items():
                    #print("occupancy:",occupancy)
                    config_indices = [config_strings.index(config_str) for config_str in configs]

                    total_density_of_occupancy_each_step=np.sum(self.crystal.ff_calculator.boundData[self.name][t_idx if len(t_idx)>1 else int(t_idx):int(t_idx+1),config_indices],axis=1) 
                    #max_total_density_of_occupancy=np.max(total_density_of_occupancy_each_step)
                    if occupancy not in self.Z() - self.charge_states:
                        # Never used.
                        self.ff_by_occupancy_and_time[self.occupancy_indices[occupancy]]=np.nan
                        self.unused_occupancy_indices.append(self.occupancy_indices[occupancy])
                        continue
                    config_ff={}
                    for config_str in configs:
                        config_idx=config_strings.index(config_str)
                        if all([self.crystal.ff_calculator.boundData[self.name][single_t_idx, config_idx]/total_density_of_occupancy_each_step[i] < 1e-3 for (i, single_t_idx) in enumerate(t_idx)]):
                            #print("Ignored",config_str)
                            #print([self.crystal.ff_calculator.boundData[self.name][single_t_idx, config_idx]/total_density_of_occupancy_each_step[i] for (i, single_t_idx) in enumerate(t_idx)])
                            config_ff[config_str]=None # don't bother calculating negligible contribution
                        else:
                            #print("PASSED",config_str)
                            #print([self.crystal.ff_calculator.boundData[self.name][single_t_idx, config_idx]/total_density_of_occupancy_each_step[i] for (i, single_t_idx) in enumerate(t_idx)])
                            #print(q_arr.shape)
                            #print(orb_occs[config_idx])
                            shell_occs=occ_dict[config_strings.index(config_str)]
                            #print(shell_occs)
                            got_ff_once=True
                            config_ff[config_str]=self.crystal.ff_calculator.ff_from_state_sane(shell_occs,q_arr,self.name)

                    all_zero_t0=None
                    for j,time in enumerate(self.times_used):
                        total_weight=0
                        ff=0
                        for config_str in configs:
                            if config_ff[config_str] is None:
                                continue
                            config_idx=config_strings.index(config_str)
                            weight=self.crystal.ff_calculator.boundData[self.name][t_idx[j], config_idx]
                            total_weight+=weight
                            ff+=weight*config_ff[config_str]
                        if total_weight==0 and j==0:
                            all_zero_t0=True and (all_zero_t0 is None or all_zero_t0) 
                            self.ff_by_occupancy_and_time[self.occupancy_indices[occupancy],j]=0
                            continue
                        if j==0:
                            all_zero_t0=False
                        assert total_weight>0, (configs, config_ff,time,occupancy,self.Z())
                        self.ff_by_occupancy_and_time[self.occupancy_indices[occupancy],j]=ff/total_weight
                        assert not np.any(np.isnan(self.ff_by_occupancy_and_time[self.occupancy_indices[occupancy],j]))
                assert got_ff_once

            #print(f"Set stochastic states for {self.name}")
        
        else:
            pass
           #self.ground_state = self.crystal.ff_calculator.get_ground_state_shells(self.name)       
    def get_num_atoms(self):
        return max(1,len(self.crystal.sym_rotations))*len(self.coords)    
    def set_coord_deviation(self,q):
        # B factor
        if self.crystal.use_bfactors:
            self.set_debye_waller(q)
        self.set_coord_deviation_no_B()
    
    def set_coord_deviation_no_B(self):
        # Random deviation for each snapshot
        num_atoms = self.num_atoms_on_coord_deviation = self.get_num_atoms()
        self.error = np.empty((num_atoms,3))
        if not self.crystal.ignore_deviations:
            for idx in range(num_atoms):
                # get random error in spherical coords based on RMS error in position, convert to cartesian.
                # there's probably a better way to do it
                err_phi,err_thet = np.random.random()*2*np.pi, np.random.random()*np.pi
                if self.name == "O" and self.crystal.random_waters is not None and idx >= self.crystal.num_non_water_oxygens:
                    err_r = np.random.normal(0, 1e3,size = (3)) # need to add this for each water atom so not same between symmetry operations
                else:
                    err_r = np.random.normal(0,self.crystal.positional_stdv, size = (3))

                self.error[idx] = err_r*(np.sin(err_phi)*np.cos(err_thet),np.sin(err_phi)*np.cos(err_thet),np.cos(err_thet))

    def get_atomic_form_factors(atom_idx,q_arr):
        '''
        (Is defined by set_atomic_form_factors).
        Returns the form factor multiplied by sqrt(I). f.shape = ( len(times) , ) + momenta.shape  
        
        atom_idx, int or int array
        q_arr = mom. transfer [1/a0], scalar or array
        '''
        raise Exception("Did not set_atomic_form_factors before calling stochastic f")
        
    def set_atomic_form_factors(self,stochastic=True):
        '''
        
        '''
        # F_i(q) = f(q)*T_i(q), where f = self.ff is the time-integrated average 
        if not stochastic:
            pass
            #self.ff = self.crystal.ff_calculator.f_average(q_arr,self.name)      # note that in taking the integral to get this ff, we included the relative intensity.
        else:
            # Undamaged case, no stochastic dynamics.
            if not self.crystal.is_damaged: 
                def get_atomic_form_factors(atom_idx,q_arr): 
                    # print("====")
                    # print("ff_times_I",self.crystal.ff_calculator.f_undamaged(q_arr,self.name)[0])
                    # print("COORDS:",self.coords)
                    # print("====")
                    return self.crystal.ff_calculator.f_undamaged(q_arr,self.name,use_intensity_for_time=self.crystal.use_intensity_for_time)[0]
            # Damaged, we 
            else:
                if self.charge_states is None:
                    def get_atomic_form_factors(atom_idx,q_arr): 
                        return self.crystal.ff_calculator.random_states_to_f_snapshots(self.times_used,self.orb_occs[atom_idx],q_arr,self.name,self.orb_occ_dict)[0]  # f has form  [times,momenta]
                else: 
                    def get_atomic_form_factors(atom_idx,q_arr): 
                        t_idx=np.searchsorted(self.crystal.ff_calculator.timeData,self.times_used)

                        occ_idx = np.zeros(shape=self.atom_occupancies[atom_idx].shape,dtype=int) - 1
                        for occ, idx in self.occupancy_indices.items():
                            occ_idx[self.atom_occupancies[atom_idx]==occ] = idx
                        if len(q_arr.shape) == 1:
                            return_val = self.ff_by_occupancy_and_time[occ_idx] * np.sqrt(self.crystal.ff_calculator.intensityData[t_idx][...,None])   # (Need to double check working as expected - not using np.vectorise)
                        elif len(q_arr.shape) == 2:
                            try:
                                return_val = self.ff_by_occupancy_and_time[occ_idx] * np.sqrt(self.crystal.ff_calculator.intensityData[t_idx][...,None,None])   # (Need to double check working as expected - not using np.vectorise)
                            except Exception as e:
                                print(self.name)
                                print(np.array(self.ff_by_occupancy_and_time).shape)
                                #print(self.atom_occupancies[atom_idx])
                                raise e
                        else:
                            assert False, f"q had shape {q_arr.shape}"
                        # print("====")
                        # print("OCC:",self.atom_occupancies[atom_idx])
                        # print("ff",self.ff_by_occupancy_and_time[self.atom_occupancies[atom_idx]])
                        # print("ff_times_I",return_val)
                        # print("COORDS:",self.coords)
                        # print("====")
                        assert not np.any(np.isnan(return_val)),(
                            np.unique(self.atom_occupancies[atom_idx]),
                            self.ff_by_occupancy_and_time[np.unique(self.atom_occupancies[atom_idx])]
                        )

                        return return_val
                       # return self.crystal.ff_calculator.avg_charge_f_snapshots(self.times_used,self.orb_occs[atom_idx],q_arr,self.name,self.orb_occ_dict,self.Z() - self.charge_states)[0]  # f has form  [times,momenta]
            self.get_atomic_form_factors = get_atomic_form_factors
class XFEL():
    def __init__(self, experiment_name, photon_energy, detector_distance_mm=100, q_minimum = None, q_cutoff = None, max_miller_idx = None, screen_type = "hemisphere", num_orients_crys=1, orientation_axis_crys = None, x_orientations = 1, y_orientations = 1, pixels_per_ring = 400, num_rings = 50,t_fineness=100,SPI_y_rotation = 0,SPI_x_rotation = 0,SPI_z_rotation = 0,all_miller_indices=False, custom_cell_dims_for_miller_indices=None,override_max_q = False,miller_indices_override=None,spot_fraction_per_orient=None):
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
            Deperecated. Determines the number of q to calculate when doing rings. Note pixels are currently just points
        t_fineness:
            Number of time steps to calculate.
        alpha,beta,gamma:
            Angle of rotation about the x, y, and z axis.
        orientation_set:
            contains each set of cardan angles for each orientation to images. Crystal only. 
            Overriden for imagings with random orientations. TODO replace x_orientations and y_orientations with this.
        x/y_orientations: 
            Number of unique x/y axis rotations to sample crystal.
        max_miller_idx:
            If specified, the maximum momentum transfer is set to that corresponding to (m,m,m), where m = max_miller_idx
        q_cutoff [1/angstrom]:
            If specified, only bragg points corresponding to momentum transfers at or below this value will be simulated.    
        """
        self.experiment_name = experiment_name
        self.detector_distance = detector_distance_mm*1e7/ang_per_bohr  # [converts to a0 (bohr)]
        self.photon_momentum = 2*np.pi/E_to_lamb(photon_energy)  # atomic units, [a0^-1]. (2pi=h)
        self.photon_energy = photon_energy  #eV Attention: not in atomic units
        self.pixels_per_ring = pixels_per_ring
        self.num_rings = num_rings
        self.num_x_orientations = x_orientations
        self.num_y_orientations = y_orientations
        self.max_miller_idx = max_miller_idx
        self.all_miller_indices = all_miller_indices
        self.miller_indices_override = miller_indices_override # overrides all selection rules TODO assert other selection rule options aren't enabled.
        self.custom_cell_dims_for_miller_indices = custom_cell_dims_for_miller_indices # This is for point of comparison with certain studies and shouldn't be used.
        self.override_max_q = override_max_q # whether should override max q when searching all miller indices 

        if self.miller_indices_override is None:
            assert spot_fraction_per_orient is None # seems unnecessary. If can't simply remove this assertion at least generate the miller indices and plug it into same logic.
        else:
            self.spot_fraction_per_orient = spot_fraction_per_orient
            if spot_fraction_per_orient is None or spot_fraction_per_orient > 1 :
                print("Defaulting to considering all Miller indices at once")
                self.spot_fraction_per_orient = 1
                
                
        if self.override_max_q:
            assert self.all_miller_indices, "Currently overriding maximum q is only supported if searching all Miller indices"
            assert self.max_miller_idx != None, "Require a maximum Miller index as maximum q is overridden." 
        if self.max_miller_idx is None:
            assert self.all_miller_indices == False, "Please specify a maximum miller index. Note that if override_max_q is false, you can choose a high value as Bragg points will still only correspond to momentum transfer below max q"
        if self.custom_cell_dims_for_miller_indices is not None:
            self.custom_cell_dims_for_miller_indices = np.array(custom_cell_dims_for_miller_indices)/ang_per_bohr
        self.min_q = 0
        if q_minimum!=None:
            self.min_q = q_minimum*ang_per_bohr

        self.hemisphere_screen = True
        eps = 0.00000000000001
        if screen_type == "flat": #well a circle
            self.hemisphere_screen = False
            self.max_q = (1-eps)*self.photon_momentum 
        elif screen_type == "hemisphere":
            # as q = 2ksin(theta), non-inclusive upper bound is 45 degrees(theoretically).  
            self.max_q = (2-eps)*self.photon_momentum/np.sqrt(2)  

        elif screen_type == "sphere":
            self.max_q = (2-eps)*self.photon_momentum
        else:
            raise Exception("Available screen types are 'flat, 'hemisphere', and 'sphere'")
        if q_cutoff != None:
            self.max_q = min(self.max_q,q_cutoff*ang_per_bohr) # convert to atomic units 

        self.t_fineness = t_fineness


        # SPI only... need to refactor
        self.phi_array = np.linspace(0,2*np.pi,self.pixels_per_ring,endpoint=False)  
        self.SPI_z_rotation = SPI_z_rotation *np.pi/180
        self.SPI_y_rotation = SPI_y_rotation *np.pi/180 # Current rotation of crystal (y axis currently) (did I mean z axis?)
        self.SPI_x_rotation = SPI_x_rotation *np.pi/180# Current rotation of crystal (y axis currently)

        ### XXX Just used for file naming purposes in scatter_MD. TODO get rid of modifying x y z rotation. 
        self.input_SPI_x_rotation=SPI_x_rotation
        self.input_SPI_y_rotation=SPI_y_rotation
        self.input_SPI_z_rotation=SPI_z_rotation
        ###

        self.z_rot_matrix = rotaxis2m(self.SPI_z_rotation,Bio_Vect(0, 0, 1))
        self.y_rot_matrix = rotaxis2m(self.SPI_y_rotation,Bio_Vect(0, 1, 0))     
        self.x_rot_matrix = rotaxis2m(self.SPI_x_rotation,Bio_Vect(1, 0, 0))   

        # crystal only...
        self.num_orientations =  num_orients_crys
        if orientation_axis_crys != None:
            if np.array(orientation_axis_crys).shape != (3,):
                raise Exception("invalid axis", orientation_axis_crys)            
            if orientation_axis_crys[0] == orientation_axis_crys[1] == orientation_axis_crys[2] == 0:
                raise Exception("axis vector has 0 length")      


        self.orientation_set = None
        if orientation_axis_crys != None:
            axis = Bio_Vect(*orientation_axis_crys)
            ori_set = [rotaxis2m(angle, axis) for angle in np.linspace(0,2*np.pi,self.num_orientations,endpoint=False)]
            self.set_orientation_set([Rotation.from_matrix(m).as_euler("xyz") for m in ori_set])
            if DEBUG or DEBUG_MODERATE:
                print("orientation set set:", self.orientation_set)
    
    def set_orientation_set(self,orientation_set):
        self.orientation_set = orientation_set
        if DEBUG or DEBUG_MODERATE:
            print("orientation set set:", self.orientation_set)
    def get_ff_calculator(self,start_time,end_time,damage_output_handle,parent_dir_path):
        ff_calculator = Plotter(damage_output_handle,parent_dir_path,out_prefix_text = "Calculating form factors...")
        plt.close()
        ff_calculator.initialise_form_factor_params(start_time,end_time,self.max_q,self.photon_energy,t_fineness=self.t_fineness) # q_fineness isn't used for our purposes.   
        return ff_calculator
    
    def fire_laser(self, start_time, end_time, sim_data_handle, sim_parent_dir_path, target : Crystal, SPI_resolution = None, results_parent_dir = RESULTS_LOCAL_PATH, circle_grid = False, pixels_across = 10, clear_output = False, random_orientation = False, SPI=False,do_not_integrate_times=False):
        """ 
        end_time: The end time of the photon capture in femtoseconds. Not a real thing experimentally, but useful for choosing 
        a level of damage. Explicitly, it is used to determine the upper time limit for the integration of the form factor.
        struct_file_path: The pdb/gromacs config file's path. Changes the variable self.atoms.
        random_orientation overrides the XFEL class's orientation_set, replacing each with a random orientation. (get same number of orientations though at present TODO.) 
        """

        print("Beginning laser")
        ff_calculator = self.get_ff_calculator(start_time,end_time,sim_data_handle,sim_parent_dir_path)
        target.set_ff_calculator(ff_calculator)    
        self.target = target
        self.integrate_times = not do_not_integrate_times

        self.used_orientations=[]


        if random_orientation == True and self.orientation_set != None:
            raise Exception("Ambiguity: random orientations set to True, but set orientations were provided.")
        if random_orientation == False and self.orientation_set == None:
            raise Exception("Providing an axis of orientation (e.g. format: [0,0,1] for z axis) or enabling random orientations is required (even with a pixel sampling simulation because I coded it bad sorry)") #TODO!
        

        if pixels_across == None and circle_grid == False and SPI:
            raise Exception("Require pixels_across argument for rectangular screen")
        
        # Create output folder for results
        directory = path.abspath(path.join(__file__ ,"../")) + "/"+ results_parent_dir + self.experiment_name + "/"  #TODO use path module properly (all should be done in function, to make end in separator have empty final arg)
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
            for rot_x in range(self.num_x_orientations):
                self.x_rot_matrix = rotaxis2m(self.SPI_x_rotation,Bio_Vect(1, 0, 0))      
                self.SPI_y_rotation = 0                 
                for rot_y in range(self.num_y_orientations):              
                    print("Imaging at x, y, rotations:",self.SPI_x_rotation,self.SPI_y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.SPI_y_rotation,Bio_Vect(0, 1, 0))      
                    self.SPI_y_rotation += 2*np.pi/self.num_y_orientations              
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
                
                    result.I += I/(self.num_y_orientations*self.num_x_orientations)
                self.SPI_x_rotation += 2*np.pi/self.num_x_orientations
                result.phi = self.phi_array
                result.q = q; assert False  # ??? this doesnt make sense. Added in assert False
            #result.package_up()
            return result
        elif SPI:
            N = pixels_across#100  NxN cells
            cell = np.zeros((N,N),dtype="object")
            result = Results_SPI()      
            
            ### Geometry
            # In crystallography, for a resolution d we have q = 2*pi/d [atomic units] as lambda=2dsin(theta), q = (4pi/lambda)sin(theta). Possibly a questionable definition for SPI without periodicity, but it is used for consistency, and Neutze 2000 makes no distinction.
                  
            # max q represents the actual limit of the change in momentum (at least in a 180 degree arc)
            # rim_q is the q we want to have at the largest unbroken ring
            rim_q = self.max_q
            if SPI_resolution!= None:
                d = SPI_resolution/ang_per_bohr # resolution
                rim_q = res_to_q(d)
                if self.max_q < rim_q:
                    print ("WARNING: resolution of " + str(SPI_resolution) + " angstroms requires q to go beyond its maximum. Using max_q="+f"{self.max_q/ang_per_bohr:.4f} (resolution {q_to_res(self.max_q/ang_per_bohr):.4f}  A) instead.")
                    rim_q = self.max_q
                print("SPI resolution:",SPI_resolution,"A", d,"a0")
            print("rim q:",rim_q*ang_per_bohr,"A", rim_q,"a0")
            #print("rim q:",rim_q)

            max_theta = self.q_to_theta(rim_q) # maximum theta of full ring.
            #screen_width = resolution_to_X(d) * 2
            screen_distance = self.detector_distance # screen-target separation [a0]
            
            # X is the distance from the centre of the screen to the point of incidence (flat screen)
            # We know where the cells are, they are equally spaced. We also know the distance from the screen to the detector.
            # This gives us theta. 
            def X_to_theta(x):
                '''
                Assumes screen distance in same units as x (a0)
                '''
                return 0.5*np.arctan2(x,screen_distance)
            def theta_to_X(theta):
                '''
                Assumes screen distance in same units as x (a0)
                '''
                return np.abs(screen_distance*np.tan(2*theta))
            
            # We must determine q at each point for finding I
            def X_to_q(x):
                '''
                Returns q [1/a0]
                '''
                theta = X_to_theta(x)
                k0 = self.photon_momentum #2*np.pi/E_to_lamb(photon_energy)
                return 2*k0*np.sin(theta)   
            
            # res. at edge corner or centre? Surely at centre, for a full ring of information
            screen_width = 2*theta_to_X(max_theta)  # edge centre  # We have placed our screen to have the desired resolution at the "rim".
            #screen_width=  2 * (np.sin(2*max_theta)/np.sqrt(2)) * screen_distance # corner
            print("Screen width:",round(screen_width*ang_per_bohr/1e7,2),"mm")            
            # Trig consistency check
            assert round(X_to_q(screen_width/2),2) == round(rim_q,2)

            ## Calculate q for each cell.
            corner_q = X_to_q(np.sqrt(2*(screen_width/2)**2))
            result.cell_width = screen_width/len(cell)
            q_grid = np.empty(cell.shape)
            result.xy = np.empty(cell.shape+(2,))
            result.I = np.zeros(cell.shape) 
            result.q = None
            result.zero_angle_I = 0 # In case dim of axis has even number of pixels           
            for i, x in enumerate(cell):
                for j, y in enumerate(x):
                    x = result.cell_width*(i-(len(cell)-1)/2)
                    y = result.cell_width*(j-(len(cell)-1)/2)
                    assert not self.hemisphere_screen, "hemisphere not supported for SPI yet"
                    
                    dist = np.sqrt(x**2+y**2)
                    q_grid[i,j] = X_to_q(dist)
                    result.xy[i,j] = np.array([x,y])
            #result.q_scr_xy = X_to_q_scr(result.xy)
            result.q_xy = X_to_q(result.xy)
            # Store a mask for the values that we should ignore.
            mask = copy.deepcopy(q_grid)
            mask[mask < self.min_q] = 0; mask[mask > self.max_q] = 0
            mask[mask != 0] = 1

            #print("MASK",mask)
            #result.I[mask == 0] = None
            result.full_ring_mask = mask
            #result.I*=mask
            # Calculate I for each cell from q (need to vectorise q)
            assert self.num_x_orientations == self.num_y_orientations==1
            for rot_x in range(self.num_x_orientations):
                self.x_rot_matrix = rotaxis2m(self.SPI_x_rotation,Bio_Vect(1, 0, 0))      
                #self.SPI_y_rotation = 0                 
                for rot_y in range(self.num_y_orientations):                  
                    print("Imaging at x, y, rotations:",self.SPI_x_rotation,self.SPI_y_rotation)
                    self.y_rot_matrix = rotaxis2m(self.SPI_y_rotation,Bio_Vect(0, 1, 0))      
                    self.SPI_y_rotation += 2*np.pi/self.num_y_orientations              
                    cell = self.generate_cell(q_grid,result.xy)
                    result.I += cell.I
                    if result.q is None:
                        result.q = cell.q
                        assert result.q.shape==result.I.shape
                    else:
                        assert result.q==cell.q
                    result.zero_angle_I += self.generate_cell(np.array([[0]]),np.array([[[0,0]]])).I #XXX

                self.SPI_x_rotation += 2*np.pi/self.num_x_orientations             # XXX (yikes)
            # convert to angstrom
            result.xy *= ang_per_bohr
            result.cell_width *= ang_per_bohr
            result.q_xy /= ang_per_bohr
            # Average out intensity # NOTE intensities aren't aligned. Shouldn't be using R factor directly on result from multiple orientations for SPI.
            # TODO Intensity should really be a list of results.
            result.I /= (self.num_y_orientations*self.num_x_orientations)


            result.save_path = directory + "SPI"+".pickle"
            result.package_up()
            return result        

        else: # Bragg reflections
            # Iterate through each orientation of crystal, picklin' up a file for each orientation
            if random_orientation:
                self.orientation_set = [(0,0,0)]*self.num_orientations  # Dummy orientations
            
            orientation_indices_override=None
            RANDOM_MILLER_FRAC = True
            if not RANDOM_MILLER_FRAC:
                if self.miller_indices_override is not None:
                    np.random.shuffle(self.miller_indices_override)
            for j, cardan_angles in enumerate(self.orientation_set):             
                
                if self.miller_indices_override is not None:
                    if self.spot_fraction_per_orient == 1:
                        orientation_indices_override = self.miller_indices_override
                    else:
                        orientation_indices_override = self.miller_indices_override.copy()
                        if RANDOM_MILLER_FRAC:
                            np.random.shuffle(orientation_indices_override)
                            orientation_indices_override = orientation_indices_override[0:int(len(self.miller_indices_override)*self.spot_fraction_per_orient)] 
                        else:
                            num_indices_total = len(self.miller_indices_override)
                            size = int(num_indices_total*self.spot_fraction_per_orient)
                            start = size*j%num_indices_total; end = start + size 
                            orientation_indices_override = orientation_indices_override[start:end] 



                print("Imaging orientation",j)
                bragg_points, miller_indices,cardan_angles = self.bragg_points(target,cell_packing = target.cell_packing, cardan_angles = cardan_angles,random_orientation=random_orientation,indices_override=orientation_indices_override)
                #print(bragg_points[3],miller_indices[3])
                #asdas

                self.used_orientations.append(cardan_angles)
                num_points = int(len(bragg_points))
                result = Results(self.target,num_points,j)
                if do_not_integrate_times:
                    result.I = np.zeros((self.t_fineness+1,)+result.I.shape) # prepend axis for time
                # Get the q vectors where non-zero
                i = 0
                #TODO vectorise
                # (Assume pixel adjacent to bragg point does not capture remnants of sinc function)
                
                max_BP =  70000 # 10920  # max number of Bragg Points to process at once 

                if len(bragg_points)>max_BP:
                    point = self.Spot(np.zeros(0),np.zeros(0),np.zeros(0))
                    num_points_left = len(bragg_points)
                    seed = np.random.randint(0,1e8) # need to seed so when iterating through each supercell the stochastic calls are the same.
                    while num_points_left > 0:
                        np.random.seed(seed)  # reset seed to same as start of orientation. TODO use generator.
                        i = len(point.q); f = len(point.q) + min(num_points_left,max_BP)
                        print(f"Iterating through Bragg points {i+1} - {f}")
                        
                        subpoint = self.generate_point(bragg_points[i:f],cardan_angles)                      
                        point.q = np.concatenate((point.q,subpoint.q),axis=0,dtype=float)
                        point.X = np.concatenate((point.X,subpoint.X),axis=0,dtype=float)
                        point.theta = np.concatenate((point.theta,subpoint.theta),axis=0,dtype=float)
                        point.phi = np.concatenate((point.phi,subpoint.phi),axis=0,dtype=float)
                        point.phi_crystal_aligned = np.concatenate((point.phi_crystal_aligned,subpoint.phi_crystal_aligned),axis=0,dtype=float)
                        point.I = np.concatenate((point.I,subpoint.I),axis=0,dtype=float)
                        num_points_left -= len(subpoint.q)
                else:
                    point = self.generate_point(bragg_points,cardan_angles)
                # fix this filling stuff
                result.phi = point.phi
                result.phi_aligned = point.phi_crystal_aligned
                result.X = point.X
                result.q = point.q
                result.I += point.I
                if do_not_integrate_times:
                    result.T = list(self.target.species_dict.values())[0].times_used # sorryyyy
                
                for_plotting = False
                if len(bragg_points)>1e4:
                    for_plotting=False
                result.package_up(miller_indices,for_plotting=for_plotting)
                #Save the result object into its own file within the output folder for the experiment
                fpath = directory + str(cardan_angles) +".pickle"
                if self.all_miller_indices:
                    version = path.basename(directory).split("_")[-1] # TODO do this in safer way
                    fpath = directory + str(j) + "_" + version + ".pickle"
                result.save_path = fpath
                with open(fpath,"wb") as pickle_out:
                    pickle.dump(result,pickle_out)

            return result
    def get_used_orientations(self):
        return self.used_orientations
    class Feature:
        def __init__(self,q,X,theta):
            self.q = q
            self.X = X # radial distance from centre of screen
            self.theta = theta          
    class Ring(Feature):
        def __init__(self,*args):
            super().__init__(*args)
    class Spot(Feature):
        def __init__(self,*args):
            super().__init__(*args)
            self.phi = np.zeros(0)
            self.phi_crystal_aligned = np.zeros(0)
            self.I = np.zeros(0)
    class Cell(Feature):
        def __init__(self,*args):
            super().__init__(*args)
            
    
    def generate_cell(self,q,xy_grid):
        solid_angle = 1  # approximate all as same for now.
        theta = self.q_to_theta(q)
        X = self.q_to_X(q)     
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")           
        cell = self.Cell(q,X,theta)
        cell.q_parr_screen = self.q_to_q_scr(q)
        phis = np.arctan2(xy_grid[:,:,1],xy_grid[:,:,0])
        #print("phis",phis)
        cell.I= self.illuminate(cell,phis=phis,SPI_proj_solid_angle = solid_angle)
        return cell

    def generate_ring(self,q,phi_array,angle_subtended):
        '''Returns the intensity(phi) array and the radius for given q.'''
        #print("q=",q)
        #r = self.q_to_r(q)# TODO
        theta = self.q_to_theta(q)
        X = self.q_to_X(q)
        solid_angle = (2*np.pi/len(phi_array)) * (np.pi/angle_subtended)   # solid angle in fractions [sp]. Solved for sphere but will be same when project to flat detector. 
        proj_solid_angle = solid_angle  # laser source.
        if self.hemisphere_screen:
            print("hemisphere/spherical screen not supported for SPI yet.")
        ring = self.Ring(q,X,theta)
        ring.I = self.illuminate(ring,phis=phi_array,SPI_proj_solid_angle=proj_solid_angle)
        return ring 
    def generate_point(self,G,cardan_angles): # G = vector
        ''' We store variables like G since they correspond to the special bragg points, unlike points from SPI that are just samples of the continuous pattern. 
        
        '''
        if len(G.shape) > 1:
            G = np.moveaxis(G,0,len(G.shape)-1)  # moves the axis corresponding to the individual momentum to the back, giving us dim = [3,num_G].
        q = np.sqrt(np.power(G[0],2)+np.power(G[1],2)+np.power(G[2],2))
        #r = self.q_to_r(q) # TODO
        X = self.q_to_X(q)
        if self.hemisphere_screen:
            X = self.q_to_q_scr_curved(G)/self.photon_momentum*self.detector_distance  # not tested...
        theta = self.q_to_theta(q)
        point = self.Spot(q,X,theta)
        point.phi = np.arctan2(G[1],G[0])
        point.G = G
        #Crystal-aligned phi. i.e. we realign all the images to be in the 0 0 0 orientation.
        G = self.rotate_G_to_orientation(G,*cardan_angles,inverse=True)[0]
        point.phi_crystal_aligned = np.arctan2(G[1],G[0])
                
        if DEBUG:
            print("X",X,"phi",point.phi,"G",G[0],G[1],G[2])


        #print("phi",point.phi,"q_parr_screen",point.q_parr_screen)

        def X_to_q(x):
            '''
            Returns q [1/a0]
            '''
            theta = X_to_theta(x)
            k0 = self.photon_momentum #2*np.pi/E_to_lamb(photon_energy)
            return 2*k0*np.sin(theta)   

        def q_to_X(self,q):
            '''
            Assumes screen distance in same units as x (a0)
            '''
            theta = self.q_to_theta(q)
            return np.abs(self.detector_distance*np.tan(2*theta))   
        def X_to_theta(x):
            '''
            Assumes screen distance in same units as x (a0)
            '''
            return 0.5*np.arctan2(x,self.detector_distance)
            



        dumb_smear_thing = False
        if dumb_smear_thing:
            xpoints = ypoints = np.linspace(0.99,1.01,3)
            G_list = []
            q_list = []
            for x in xpoints:
                for y in ypoints:
                    G_copy = G.copy()
                    G_copy[0]*=x
                    G_copy[1]*=y
                    G_copy[2]*=x*y  # this is so dumb ugh 
                    G_list.append(G_copy)
                    q_list.append(np.sqrt(np.power(G_copy[0],2)+np.power(G_copy[1],2)+np.power(G_copy[2],2)))  
            point.I = None
            for _G, _q, in zip(G_list,q_list):
                pointCopy = self.Spot(_q,X,self.q_to_theta(_q))
                pointCopy.phi = np.arctan2(_G[1],_G[0])
                pointCopy.G = _G
                if point.I is None:
                    point.I = self.illuminate(pointCopy,cardan_angles=None)
                else:
                    point.I += self.illuminate(pointCopy,cardan_angles=None)
        else:
            point.I = self.illuminate(point,cardan_angles=cardan_angles)




        #Using formula in caleman 2011 but arbitrary scale
        #lamb =  E_to_lamb(self.photon_energy)
        #lorentzFactor = lamb**3/np.sin(point.phi)
        #lorentzCorrection = 1/np.sin(2*self.q_to_theta(point.q))
        #point.I*=lorentzCorrection
        
        return point   
    




    # Returns the relative intensity at point q for the target's unit cell, i.e. ignoring crystalline effects.
    # If the feature is a bragg spot, this gives its relative intensity, but due to photon conservation won't be the same as the intensity without crystallinity - additionally different factors for non-zero form factors occur across different crystal patterns.
    def illuminate(self,feature, phis = None,cardan_angles = None, SPI_proj_solid_angle=None,seed=None):  # Feature = ring or spot.
        """Returns the intensity at q. Not crystalline yet.
        phis is an ndarray containing the angle of each point to calculate for (SPI) rings, but should contain only 1 element for points.
        """
        # if phis == None:
        #     phis = self.phi_array


        if seed is not None:
            np.random.seed(seed)
        
        if type(feature) == self.Spot:
            SPI = False
        elif type(feature) in [self.Cell,self.Ring]:
            SPI = True

        F_shape = tuple()
        if type(feature) is self.Spot:
            F_shape += (self.t_fineness+1,)           
            if len(feature.G.shape) > 1:
                F_shape += (len(feature.G[2]),) #[times,num_G]   
        else:
            if type(feature) is self.Ring:
                F_shape += phis.shape  
            F_shape += (self.t_fineness+1,)# [?phis?,times]
            if type(feature.q) is np.ndarray:
                F_shape += feature.q.shape          # [?phis?,times,feature.q.shape]
        F_supercells = np.zeros(self.target.supercell_simulations,dtype="object")
        non_empty_species_dict:dict[str,Atomic_Species] = {}
        for k, v in self.target.species_dict.items():
            if len(v.coords) >0:
                non_empty_species_dict[k]=v
        for species in non_empty_species_dict.values():
            species.ff_by_occupancy_and_time=None
        for S in range(self.target.supercell_simulations):   
            print("Simulating supercell", S)
            times_used = None
            self.target.set_stochastic_positions(feature.q)
            for species in non_empty_species_dict.values():
                species.set_stochastic_electronic_states(q_arr=feature.q)   
                species.set_atomic_form_factors()
                if times_used is None:
                    times_used = species.times_used
                else:
                    if DEBUG:
                        assert np.all(times_used == species.times_used)
            if (times_used[-1] == times_used[0] and times_used.size!=1):   
                raise Exception("Intensity array's final time equals its initial time")                  
            # Technically sum of F(t)*sqrt(J(t)), where F = sum(f(q,t)*T(q)), and J(t) is the incident intensity, thus accounting for the pulse profile. (J(t) is accounted for in get_atomic_form_factors)
            F_sum = np.zeros(F_shape,dtype="complex_")  
            for species in non_empty_species_dict.values():
                print(f"Atom {species.name}")
                if DEBUG or DEBUG_MODERATE:
                    print("------------------------------------------------------------")
                    print("Getting contribution to integrand from species",species.name)
                if not np.array_equal(times_used,species.times_used):
                    raise Exception("Times used don't match between species.")        
                # iterate through every atom including in each symmetry of unit cell (each asymmetric unit)
                max_atoms_per_loop = 20 # Restrict array size to prevent computer explosions. 
                self.target.reset_unique_points()
                for s in range(len(self.target.sym_rotations)):
                    #print("Working through symmetry",s)
                    num_atom_batches = int(len(species.coords)/max_atoms_per_loop)+1
                    global inner_loop
                    def inner_loop(a_batch):
                        F_sum = np.zeros(F_shape,dtype="complex_")
                        if a_batch%100==0:
                            print(f"atom batch {a_batch}/{num_atom_batches}")
                        atm_idx = np.arange(len(species.coords)*s+max_atoms_per_loop*a_batch, len(species.coords)*s + min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                        if len(atm_idx) == 0: 
                            return F_sum
                        relative_atm_idx = np.arange(max_atoms_per_loop*a_batch, min(len(species.coords), max_atoms_per_loop*(a_batch+1)))
                        R = np.array(species.coords[relative_atm_idx[0]:relative_atm_idx[-1]+1]) 
                        coord = self.target.get_sym_xfmed_point(R,s)
                        if not self.target.ignore_deviations:# or self.target.use_bfactors:
                            coord += species.error[relative_atm_idx] # dim = [ N, 3], where N is number of coords.
                        if DELETENONUNIQUE:
                            unique_mask = self.target.get_sym_xfmed_point_unique_mask(R,s)
                            coord = coord[unique_mask]
                            atm_idx = atm_idx[unique_mask]
                            relative_atm_idx = relative_atm_idx[unique_mask]
                            #print(species.name,coord.shape)

                        if SPI:
                            # Rotate to target's current orientation 
                            # rot matrices are from bio python and are LEFT multiplying. TODO should be consistent replace this with right mult. 
                            coord = coord @ self.x_rot_matrix  
                            coord = coord @ self.y_rot_matrix
                            coord = coord @ self.z_rot_matrix
                        # Get spatial factor T
                        if SPI:
                            T = self.SPI_interference_factor(phis,coord,feature)  # if grid: [phis,qx,qy]  if ring: [phis,q] (unimplemented)
                        else:
                            T= self.interference_factor(coord,feature,cardan_angles)  #[num_G] 
                        if self.target.use_bfactors and not self.target.zero_bfactors:
                            T*=species.debye_waller_factor[relative_atm_idx]
                        f = species.get_atomic_form_factors(atm_idx, feature.q)  / np.sqrt(self.target.num_cells*self.target.num_supercells) # Dividing by np.sqrt(self.num_cells) so that fluence is same regardless of num cells. 
                        if self.miller_indices_override is None:
                            # resolution limit imposed by wavelength 
                            f*=(feature.q < self.max_q)
                        assert np.sqrt(self.target.num_cells*self.target.num_supercells)>0
                        assert not np.any(np.isnan(f)),f
                        # print("##########")
                        # print("f",f)
                        # print("T",T)
                        # print(len(self.target.sym_rotations))
                        # print("##########")
                        
                        #print(F_sum.shape,T.shape,f.shape)                             
                        if SPI: 
                            if type(feature) is self.Cell:
                                F_sum += np.sum(T[:,None,...]*f,axis=0)          #[num_atoms,None,qX,qY] X [num_atoms,times,qX,qY] -> [times,qX,qY]
                            if type(feature) is self.Ring:           
                                F_sum += np.sum(T[:,:,None,:] * f[:,None,:,:],axis=0)        # [num_atoms,phis,None,num_|q|]X[num_atoms,None,times, num_|q|]  -> [phis,times,num_|q|]
                        else:
                            F_sum += np.sum(T[:,None,:] * f,axis=0)                           # [num_atoms,None,num_G]X[num_atoms,times,num_G]  ->[times,num_G] 
                        #print("s,F_sum",s,F_sum)
                                                                          #I =  np.square(np.abs(F_sum[:,0]))  # 
                        return F_sum
                    with Pool(NUM_THREADS,maxtasksperchild=100) as p:
                        # Important that we iterate, and don't convert to list (e.g. F_sum+=np.sum(list(F_sums),axis=0)) as would be very big memory allocation
                        for other_Fsum in p.map(inner_loop,range(num_atom_batches)):
                            F_sum+=other_Fsum  
                        
            F_supercells[S] = F_sum
        # Add supercells
        F_cry = np.zeros(F_shape,dtype="complex_")

        length = np.ceil(self.target.num_supercells**(1/3)) 

        x, y, z= np.meshgrid(np.arange(0, length), np.arange(0, length), np.arange(0, length))
        super_cube_coords = np.stack([ y.flatten(), x.flatten(), z.flatten()], axis = -1) # sadly not dimensionally-transcendental enough to be prefixed "hyper"        
        curr_cube_idx = 0
        super_batch_size = 50
        supercells_remaining = self.target.num_supercells
            
        if self.target.cell_packing == "triclinic":
            triclinic_basis = get_triclinic_basis(self.target.cell_angles)
        while supercells_remaining > 0:
            end_cube_idx = curr_cube_idx + min(supercells_remaining,super_batch_size)
            super_coords = np.empty((self.target.num_supercells,3))  
            super_coords = super_cube_coords[curr_cube_idx:end_cube_idx+1]*self.target.supercell_dim
            assert (super_coords.shape[0]>0), f"{super_coords.shape[0]} {super_coords.shape} {curr_cube_idx} {end_cube_idx}"
            assert(super_coords.shape[0]>0)
            if self.target.cell_packing == "triclinic":
                super_coords= (super_coords*self.target.supercell_dim)@triclinic_basis


            F_supercell_copies = np.zeros(shape=(len(super_coords),)+(F_supercells[0].shape),dtype=complex)
            for p in range(len(F_supercell_copies)):
                F_supercell_copies[p] = np.random.choice(F_supercells)
            if self.all_miller_indices:
                # Assume at Bragg conditions... I THINK THIS IS WRONG WE WILL NEED TO ROTATE G TO BRAGG CONDITION FOR EACH AND AVERAGE OVER.
                T_supercell = np.ones((super_coords.shape[0],feature.G.shape[-1])) # Only constructive interference at bragg spots.
            else:
                if SPI:
                    T_supercell = self.SPI_interference_factor(phis,super_coords,feature)
                else:
                    T_supercell = self.interference_factor(super_coords,feature,cardan_angles,rotate_back=False)
            F_cry += np.sum(F_supercell_copies*T_supercell[:,None],axis=0)
            supercells_remaining -= super_batch_size
            curr_cube_idx = end_cube_idx
        
        I = np.square(np.abs(F_cry))
        I_ref, _ = self.target.ff_calculator.I_avg()
        I/=1e15 
        I_ref/=1e15
        if self.integrate_times:
            if times_used.size>1:
                # Integrate over time to get the intensity
                time_axis = 0
                if type(feature) is self.Ring:
                    time_axis = 1
                I = np.trapz(I,times_used,axis = time_axis) / (times_used[-1]-times_used[0])       #[num_G] for points, or for SPI: [phis,feature.q.shape], corresponding to rings or square grid
            else:
                I = I[0]  

        # For unpolarised light (as used by Neutze 2000). (1/2)r_e^2(1+cos^2(theta)) is thomson scattering - recovering the correct equation for a lone electron, where |f|^2 = 1 by definition.    
        # Generally not important due to small angles involved.
        r_e_sqr = 2.83570628e-9
        I*= r_e_sqr*(1/2)
        I_ref*= r_e_sqr*(1/2)
        if not self.all_miller_indices:  # TODO maybe turn this off in gneeral since we are assuming it is corrected for
            thet = self.q_to_theta(feature.q)
            I*=(1+np.square(np.cos(2*thet))) 
        else:
            I/=I.size # divvied up by spots

        assert not np.any(np.isnan(I)),I

        if SPI:
            # For crystal we can approximate infinitely small pixels and just consider bragg points.
            # But for SPI need to take the pixel's size into account. Neutze 2000 makes the following approximation:
            I *=  SPI_proj_solid_angle    # equiv. to *= solid_angle
        
        assert not np.any(np.isnan(I)), (I,SPI_proj_solid_angle)

        print("Total screen-incident intensity TODO not matching below with all miller  = ","{:e}".format(np.sum(I)))
        print("Intensity scattered by free electron TODO not matching above with all miller = ","{:e}".format(I_ref))
        return I


    def interference_factor(self,coord,feature,cardan_angles,rotate_back=True): # TODO Remove cardan_angles input
        """ theta = scattering angle relative to z-y plane """ 
        if DEBUG:
            assert coord.shape[0]>0
        q_vect = feature.G.copy()
        #Rotate our G vector BACK to the real laser orientation relative to the crystal.
        if rotate_back: # Should still be same result regardless of orientation TODO double check
            q_vect = self.rotate_G_to_orientation(feature.G.copy(),*cardan_angles,inverse=True)[0]     
        coord = np.moveaxis(coord,0,-1)  #  dim = [xyz,atoms]
        q_vect = np.moveaxis(q_vect,0,-1) # dim = [momenta,xyz]
        #q_dot_r = np.apply_along_axis(np.dot,len(q_vect.shape)-1,q_vect,coord) # dim = [num_G]  
        q_dot_r = np.apply_along_axis(np.dot,len(q_vect.shape)-1,q_vect,coord) # dim = [num_G]  

        q_dot_r = np.moveaxis(q_dot_r,-1,0)
        coord = np.moveaxis(coord,-1,0)     
        T = np.exp(-1j*q_dot_r)
        # T = np.exp(0*q_dot_r)
        # print(T)
        # adsads
        return T  # [Num_atoms, num_G (bragg spots)]
    
    def random_water_time_interference_factor(self,coord,feature,cardan_angles,times,target):
        q_vect = self.rotate_G_to_orientation(feature.G.copy(),*cardan_angles,inverse=True)[0]
        T = np.zeros((coord.shape[0],len(times),q_vect.shape[1]))
        target.random_water_start_idx
        for t in times:
            target.reinitialize_random_waters()
        assert False, "TODO" #TODO

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
    def bragg_points(self,crystal, cell_packing, cardan_angles,random_orientation=False,indices_override=None):
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
        # a = primitive lattice vectors, b = reciprocal lattice vectors
        if cell_packing == "SC" or cell_packing == "triclinic" or cell_packing == "FCC" or cell_packing == "BCC" or cell_packing == "FCC-D":
            if cell_packing == "SC":
                a = np.array([[1,0,0],[0,1,0],[0,0,1]])
            if cell_packing == "BCC":
                a = 0.5*np.array([[-1,1,1],[1,-1,1],[1,1,-1]])
            if cell_packing == "FCC":
                a = 0.5*np.array([[0,1,1],[1,0,1],[1,1,0]])
            if cell_packing == "triclinic":
                a =get_triclinic_basis(crystal.cell_angles)
                
            if self.custom_cell_dims_for_miller_indices is None:
                a = np.multiply(a.T,crystal.cell_dim).T
            else:
                print(f"WARNING, using custom cell dims {self.custom_cell_dims_for_miller_indices}")
                a = np.multiply(a,self.custom_cell_dims_for_miller_indices) 
        else: 
            raise Exception("Unknown cell packing type")
        if not np.array_equal(crystal.supercell_dim,crystal.cell_dim):
            print("Warning: This code samples the intensity at the miller indices, and thus does not capture peak broadening - consider using the SPI imaging mode instead.")
        b1 = np.cross(a[1],a[2])
        b2 = np.cross(a[2],a[0])
        b3 = np.cross(a[0],a[1])
        b = 2*np.pi*np.array([b1,b2,b3])/(np.dot(a[0],np.cross(a[1],a[2])))
        
        assert not np.isnan(b).any() , f"primitive vectors {a} reciprocal vectors {b}"
        
        ###################### NOT TESTED TODO
        if indices_override is not None:
            G, cardan_angles = get_G(indices_override); 
            indices = indices_override
        else:
        ######################

            # Cast a wide net, catching all possible permutations of miller indices.
            #   G = hb1 + kb2 + lb3. (bi = lattice vector)
            #   !Attention! Assuming vectors are orthogonal.
            # TODO double check not cutting off possible values.
            if self.override_max_q:
                h_max = k_max = l_max = self.max_miller_idx    
            else:
                q1 = sum(pow(self.max_miller_idx*element, 2) for element in b[0])
                q2 = sum(pow(self.max_miller_idx*element, 2) for element in b[1])
                q3 = sum(pow(self.max_miller_idx*element, 2) for element in b[2])
        
                print(res_to_q(self.max_q/ang_per_bohr))
                if DEBUG:
                    if self.max_miller_idx:
                        highest_possible_q = np.sqrt(q1+q2+q3)/ang_per_bohr
                        print(f"q for ({[self.max_miller_idx,]*3}): {highest_possible_q}, resolution: {q_to_res(highest_possible_q)}")
                    
                q_1_max = self.max_q
                q_2_max = self.max_q
                q_3_max = self.max_q        # q = (0,0,l) case.
                h_max = 0
                while (np.sqrt(sum(pow(h_max*element, 2) for element in b[0])) <= q_1_max 
                    and h_max < self.max_miller_idx):
                    h_max += 1
                k_max = 0
                while (np.sqrt(sum(pow(k_max*element, 2) for element in b[1])) <= q_2_max
                    and k_max < self.max_miller_idx):
                    k_max += 1 
                l_max = 0
                while (np.sqrt(sum(pow(l_max*element, 2) for element in b[2])) <= q_3_max
                    and l_max < self.max_miller_idx):
                    l_max += 1     

                q1 = sum(pow(h_max*element, 2) for element in b[0])
                q2 = sum(pow(k_max*element, 2) for element in b[1])
                q3 = sum(pow(l_max*element, 2) for element in b[2])
                highest_possible_q = np.sqrt(q1+q2+q3)
                if DEBUG or DEBUG_MODERATE:
                    print(f"q for ({h_max},{k_max},{l_max}): {highest_possible_q/ang_per_bohr}, resolution: {q_to_res(highest_possible_q/ang_per_bohr)}")
                assert(highest_possible_q > self.min_q), f"highest allowed q{ {highest_possible_q/ang_per_bohr}} < smallest q {self.min_q/ang_per_bohr}"

            h_set = np.arange(-h_max,h_max+1,1)
            k_set = np.arange(-k_max,k_max+1,1)
            l_set = np.arange(-l_max,l_max+1,1)
            indices = list(itertools.product(h_set,k_set,l_set))
            indices = np.array([*set(indices)])   
        
            # Selection rules
            if cell_packing == "SC" or cell_packing == "triclinic":
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

            def miller_selection_rule():
                return None
            if self.max_miller_idx != None:
                m = self.max_miller_idx
                max_g_vect = get_G(np.full((1,3),m))[0][0]
                if not self.override_max_q:
                    self.max_q = min(self.max_q,np.sqrt(((max_g_vect[0])**2+(max_g_vect[1])**2+(max_g_vect[2])**2)))
                else:
                    self.max_q = np.sqrt(((max_g_vect[0])**2+(max_g_vect[1])**2+(max_g_vect[2])**2))
                def miller_selection_rule(indices): # Probably unnecessary I'm just making sure...
                    return  (abs(indices[0]) <= m and abs(indices[1]) <= m and abs(indices[2]) <= m)
            #print("max q (i.e. rim q):",self.max_q/ang_per_bohr)
            
            print("using q range of ", self.min_q/ang_per_bohr,"-",self.max_q/ang_per_bohr," angstrom-1")
            print("corresponding to max resolution: ", q_to_res(self.max_q)*ang_per_bohr," angstrom")
            def min_max_q_rule(g):
                return self.min_q <= np.sqrt(((g[0])**2+(g[1])**2+(g[2])**2)) <= self.max_q
            #max_q_rule = lambda f: np.sqrt(((f[0]*np.average(cell_dim))**2+(f[1]*np.average(cell_dim))**2+(f[2]*np.average(cell_dim))**2))<= self.max_q
            def select_friedel(indices):
                return np.all(indices>=0)
            # Catch the miller indices with a boolean mask
            if self.all_miller_indices:
                mask= np.apply_along_axis(miller_selection_rule,1,indices)
                mask *= np.apply_along_axis(selection_rule,1,indices)
                mask *= np.apply_along_axis(select_friedel,1,indices)
                #TODO select for just one symmetry.
                if not self.override_max_q:
                    mask*=np.apply_along_axis(min_max_q_rule,1,G_temp)
            else:
                mask = np.apply_along_axis(selection_rule,1,indices)*np.apply_along_axis(min_max_q_rule,1,G_temp)*np.apply_along_axis(self.mosaic_elastic_condition,1,G_temp)
                if self.max_miller_idx != None:
                    mask*=np.apply_along_axis(miller_selection_rule,1,indices)
            indices = indices[mask]
        assert(len(indices) > 0)


        actual_max_q = 0
        max_q_indices= [0,0,0]
        for h,k,l in indices:
            q1 = sum(pow(h*element, 2) for element in b[0])
            q2 = sum(pow(k*element, 2) for element in b[1])
            q3 = sum(pow(l*element, 2) for element in b[2])
            q = np.sqrt(q1+q2+q3)
            if actual_max_q < q:
                actual_max_q = q
                max_q_indices = h,k,l
        assert max_q_indices != [0,0,0]
        print(f"Best resolution point {max_q_indices}: {q_to_res(actual_max_q)*ang_per_bohr} angstrom." )



        print("Cardan angles:",cardan_angles)
        print("Number of points:", len(indices))   
        # Commented out because now break up iteration indices if too large.
        # assert len(indices) < 10920, "Output size failsafe triggered, output size likely >~ 1 GiB."
        # if len(indices) > 2000:
        #     print("WARNING: very high number of points!")        
        for elem in indices:
            if DEBUG:
                print(elem)
        

        G, cardan_angles = get_G(indices)
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
            """ Returns the screen-parallel component of q - useful for getting q_x and q_y
            """
            theta = self.q_to_theta(q)
            return q/(np.sin(np.pi/2-theta))      
    def q_to_X(self,q):
        '''
        Assumes screen distance in same units as x (a0)
        '''
        theta = self.q_to_theta(q)
        return np.abs(self.detector_distance*np.tan(2*theta))    
               
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

# Rotate to non-orthogonal axis
def get_triclinic_basis(primitive_angles): # angles in degrees
    # is this the transverse?
    #a =  np.array([[1,0,0],[0,1,0],[0,0,1]]).astype(np.float64)
    A,B,C = np.deg2rad(primitive_angles)

    c_1 =  np.cos(B)
    c_2 = (np.cos(A)-np.cos(B)*np.cos(C))/np.sin(C)
    c_3 = np.sqrt(1-c_1**2-c_2**2)

    return np.array(
        [[1,0,0],
        [np.cos(C),np.sin(C),0],
        [c_1,c_2,c_3]]
    )
    #for i,  angle in zip(range(len(a)),primitive_angles):
        #angles = np.array([0,0,0],dtype=np.float64)
        #angles[i] = angle-90
        #angles = np.array(primitive_angles) - 90
        
        #a_vector = np.array([0,0,0])
        #a_vector = a[i]
        #a[i] = #a_vector @ Rotation.from_euler('xyz',angles,degrees=True).as_matrix()

        #a[0] = np.cos() + np.sin()
    #return a

def E_to_lamb(photon_energy):
    """Energy (eV) to wavelength in A.U."""
    E = photon_energy  # eV
    return 2*np.pi*c_au/(E/eV_per_Ha)
        #
        # q = 4*pi*sin(theta)/lambda = 2pi*u, where q is the momentum in AU (a_0^-1), u is the spatial frequency.
        # Bragg's law: n*lambda = 2*d*sin(theta). d = gap between atoms n layers apart i.e. a measure of theoretical resolution.

def scatter_scatter_plot(get_R_only = False,neutze_R = True, crystal_aligned_frame = False ,SPI_result1 = None, SPI_result2 = None, full_range = True,num_arcs = 50,num_subdivisions = 40, result_handle = None, results_parent_dir = RESULTS_LOCAL_PATH, compare_handle = None, normalise_intensity_map = False, show_grid = False, cmap_power = 1, cmap = None, min_alpha = 0.05, max_alpha = 1, 
                         bg_colour = "grey",solid_colour = "white", show_labels = False, radial_lim = None, plot_against_q=False,log_I = True, log_dot = False,  fixed_dot_size = False, dot_size = 1, crystal_pattern_only = False, log_radial=False,cutoff_log_intensity = None,
                         spi_full_rings_only=True,min_R_dmg_pixel = 0.1,cmap2=None,log_range=None,custom_fig_width=None,custom_fig_height=None,log_diff_vmin = -0.5, log_diff_vmax = 0.5, normalize_to_centre = True, dpi=100,
                         plot_handle="",show_plot=True):
    ''' (Complete spaghetti at this point.)
    Plots the simulated scattering image.
    result_handle:

    compare_handle:
        results/compare_handle/ is the directory of results that will be subtracted from those in the results directory.
        Must have same orientations.
    
    '''
    # Pixel plots...
    if custom_fig_height is None:
        custom_fig_height = plt.rcParams['figure.figsize'][1]
    if custom_fig_width is None:
        custom_fig_width = plt.rcParams['figure.figsize'][0]
    plt.rcParams['figure.figsize'] = [custom_fig_width,custom_fig_height]
    plt.rcParams['figure.dpi'] = 800
    LOG10 = False
    if LOG10:
        log_function = np.log10
        if log_range is None:
            log_range = 5
    else:
        log_function = np.log
        if log_range is None:
            log_range = 10
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
                bad_val = False
                mod_val = val
                if compare_handle == None:
                    # Plot actual image!
                    for elem in (val,min_z,max_z,cmap_power):
                        if elem == -np.inf or elem == None or elem == np.inf:
                            bad_val = True   
                    if val < min_z or val > max_z:
                        bad_val = True
                    if bad_val:
                        raise Exception("invalid value within",val,min_z,max_z,cmap_power)
                    mod_val = ((val-min_z)/(max_z-min_z))**cmap_power
                
                if mod_val < 0 or mod_val > 1:
                    raise Exception("rgb val outside acceptable range")
                return self.scalarMap.to_rgba(mod_val)            
    else:
        # Broken TODO
        r,g,b = to_rgb(solid_colour)
        colours = [(r,g,b,a) for a in np.clip(z/max_z,min_alpha,max_alpha)]
    if cmap2 is None:
        cmap2 = cmap
    #def add_screen_properties(fig_width=14,fig_height=8.4):
    def add_screen_properties(fig_width=None,fig_height=None):
        if fig_width is None:
            fig_width = custom_fig_width
        if fig_height is None:
            fig_height = custom_fig_height
        if radial_lim:
            bottom,top = plt.ylim()
            plt.ylim(bottom,radial_lim)
        if log_radial:
            plt.yscale("log")  
        plt.gca().set_facecolor(bg_colour)      
        #plt.gcf().set_figwidth(20)         # if not using widget magic
        #plt.gcf().set_figheight(20)        #       
            
        plt.gcf().set_figwidth(fig_width)
        plt.gcf().set_figheight(fig_height)                  


    # (If Bragg spots)
    if result_handle != None: #TODO replace this atrocious way of distinguishing between inf. crystal and finite
        if plot_against_q:
            radial_lim /= ang_per_bohr
        else:
            radial_lim*= ang_per_bohr        
        
        results_dir = results_parent_dir+result_handle+"/"
        compare_dir = None
        if compare_handle!= None:
            compare_dir = results_parent_dir+compare_handle+"/"        

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


        def plot_spots(I_ideal,I_real, miller_indices,normalise=False):
            plt.close()
            stringified = []
            for elem in miller_indices:
                stringified.append(str(elem[0])+str(elem[1])+str(elem[2])) 
            if normalise:
                I_real *= np.sum(I_ideal)/np.sum(I_real) #normalise
            plt.bar(stringified,np.sqrt(I_ideal),alpha=1)
            plt.bar(stringified,np.sqrt(I_real),alpha=1,color='r',width=0.4)
            plt.ylim(0,np.sqrt(max(np.max(I_ideal),np.max(I_real))))
            plt.xticks(rotation="vertical")
            plt.show()
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
            max_z = log_function(max_z)
            min_z = log_function(min_z)
        if DEBUG or DEBUG_MODERATE:
            if log_I:
                print("log(I) max,min ",max_z,min_z)
            else:
                print("I max,min",max_z,min_z)
        if max_z == min_z:
            print("Single intensity detected. Ignoring max/min z")
            min_z = max_z-1
        #TODO dont even use min z?

        # Plot each orientation's scattering pattern/add each to sector histogram
        added_colorbar = False
        for_plotting = True
        for filename in os.listdir(results_dir):
            result1,result2 = get_result(filename,results_dir,compare_dir)
            if result1 == "__PASS__":
                continue
            if result1 == None:
                break         
            # get points at same position on screen.
            radial_axis = result1.X*ang_per_bohr/1e7
            if plot_against_q:
                radial_axis = result1.q_scr/ang_per_bohr
            radial_axis = radial_axis[0]    

            # Plot spot intensities
            plot_all_the_spots = True
            if plot_all_the_spots:
                if result2 != None and crystal_aligned_frame:
                    plot_spots(result2.I.flatten(), result1.I.flatten(),result1.miller_indices)


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
            if result1.for_plotting:
                for i in range(len(result1.I)):
                    if (radial_axis[i], phi[i],result1.image_index)  in processed_copies:
                        unique_values_mask[i] = False
                        continue
                    else:
                        unique_values_mask[i] = True
                    # if (radial_axis[i] > radial_lim):
                    #     unique_values_mask[i] = False
                        
                    # Average out spots (TODO need to check they are the same for no stochastic variation.)
                    # TODO need to change this to be the bragg indices.
                    matching_rad_idx = np.nonzero(radial_axis == radial_axis[i])  # numpy note: equiv. to np.where(condition). Non-zero part irrelevant.
                    matching_phi_idx = np.nonzero(phi == phi[i])
                    for elem in matching_rad_idx[0]:
                        if elem in matching_phi_idx[0]:
                            identical_count[i] += 1
                            matching1 = result1.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                            for value1 in matching1:
                                tmp_I1[i] += value1    
                            if result2 != None:
                                matching2 = result2.I[(radial_axis == radial_axis[i])*(phi == phi[i])]
                                for value2 in matching2:
                                    tmp_I2[i] += value2      
                #print(identical_count)                                                
                tmp_I1 /= (identical_count)
                tmp_I2 /= (identical_count)

                    #processed_copies.append((radial_axis[i],phi[i],result.image_index[i]))
                #print("processed copies:",processed_copies)

                # Get the dependent variable.
                result = copy.deepcopy(result1)
                result.I = tmp_I1 
                # get z used for colour of scatter plot.
                if compare_dir != None:
                    result2.I = tmp_I2                    
                    result.diff(result2)
                    z = result.R[unique_values_mask]   # Making comparison, set z to be measure of difference 
                else:
                    z = result.I[unique_values_mask]   # no comparison, z is intensity

                identical_count = identical_count[unique_values_mask]
                # Intensities used for histogram 
                I1 = tmp_I1[unique_values_mask]
                I2 = tmp_I2[unique_values_mask]
                radial_axis = radial_axis[unique_values_mask]
                phi = phi[unique_values_mask] + np.pi/2 # to align with the SPI pixel plot
                phi[phi>=np.pi] -= 2*np.pi
            else:
                results_for_plotting = False

            #sector_histogram += get_histogram_contribution(z,result.phi,radial_axis)
            if result2 != None and results_for_plotting:
                ## Sectors/Fragmented rings ##
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis)
                # Add to sum of all orientations' results
                sector_num_histogram += numerators
                sector_den_histogram += denominators                

                ## Full rings ## Note: we are filling 2d arrays with same-valued arcs for each radius so that pcolormesh can plot circles. 
                # The numerators/denominators are effectively tiled, e.g.: numerators = np.tile(numerators,(1,num_arcs))
                numerators,denominators = get_sector_histogram_contribution(I1,I2,phi,radial_axis,phi_edges = np.array([-np.pi,np.pi]))
                R_num_histogram += numerators
                R_den_histogram += denominators            
                # if neutze_R:
                    
                #     all_I_real.extend(I1)
                #     all_I_ideal.extend(I2)

                    


    

            #debug_mask = (0.01 < radial_axis[0])*(radial_axis[0] < 100)
            #print(identical_count[debug_mask])
            #print(radial_axis[0][debug_mask])
            #print(phi[debug_mask ]*180/np.pi)
            if not get_R_only and results_for_plotting:
                
                if log_I: 
                    z = log_function(z)  
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
                alpha_modified_cmap_colours = np.empty((len(z),4)) 
                for i, K in enumerate(z):
                    try:
                        rgba = COL.get_rgb(K)
                    except Exception as e:
                        raise Exception("(max/min z:" +str(max_z) + "/" + str(min_z) + ") - val " + str(K) + " did not work for get_rgb().Original error: " + str(e))     
                    #Replace alpha
                    rgba=rgba[0:3] + (np.clip((K*(max_alpha-min_alpha))/(max_z) + min_alpha,min_alpha,None),)
                    alpha_modified_cmap_colours[i] = np.array(rgba)                           
                

                colours = [(r,g,b,a) for r,g,b,a in alpha_modified_cmap_colours] 
                colours = np.around(colours,10)   

                if not added_colorbar:
                    if compare_dir == None: 
                        if DEBUG or DEBUG_MODERATE:
                            print("Warning, colorbar not working with color power at present. Need to create cmap from COL")
                        if not normalise_intensity_map:
                            # Good for debugging
                            fig.colorbar(cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=np.min(z),vmax=np.max(z)),cmap=cmap),ax=ax)
                            if DEBUG or DEBUG_MODERATE:
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

                sc = ax.scatter(phi,radial_axis,c=colours,s=s)
                plt.grid(alpha = min(show_grid,0.6),dashes=(5,10))   
                if show_labels:
                    rad_max = 0
                    q_max_idx = None
                    max_miller_index = 0
                    q_max = 0
                    for i, miller_indices in enumerate(result.miller_indices[unique_values_mask]):
                        #if rad_max < radial_axis[i]:
                        if q_max < result1.q[0][unique_values_mask][i]/ang_per_bohr:
                            rad_max = radial_axis[i]
                            q_max = result1.q[0][unique_values_mask][i]/ang_per_bohr
                            q_max_idx = i
                        max_miller_index = max(max_miller_index,np.max(miller_indices))
                        ax.annotate("%.0f" % miller_indices[0]+","+"%.0f" % miller_indices[1]+","+"%.0f" % miller_indices[2], (phi[i], radial_axis[i]),ha='center')
                    print("Num points", len(result.miller_indices[unique_values_mask]))
                    print("Max miller index",max_miller_index)
                    print("Miller indices of max q =","%.0f" % result.miller_indices[unique_values_mask][q_max_idx][0]+","+"%.0f" % result.miller_indices[unique_values_mask][q_max_idx][1]+","+"%.0f" % result.miller_indices[unique_values_mask][q_max_idx][2]) 
                    if plot_against_q:
                        print("max q_scr =",rad_max)
                    else:
                        print("max radius=",rad_max)
                    print("q max=", q_max)

        # Get merged intensities with miller indices in each set
        miller1, I1_all = read_scalepack(result_handle)
        miller2, I2_all = read_scalepack(compare_handle)
        miller_indices_1  = np.empty(shape=miller1.shape[0],dtype=object)
        miller_indices_2  = np.empty(shape=miller2.shape[0],dtype=object)

        for i, m in enumerate(miller1):
            miller_indices_1[i] = list(m)
        for i, m in enumerate(miller2):
            miller_indices_2[i] = list(m)
        

        miller_indices = np.intersect1d(miller_indices_1,miller_indices_2)
        print(f"{miller_indices.shape[0]} Bragg spots compared")
        I1 = np.zeros(miller_indices.shape)
        I2 = np.zeros(miller_indices.shape)
        for i, hkl in enumerate(miller_indices):
            hkl = list(hkl)
            if DEBUG:
                print(f"1 num points for {hkl}",np.count_nonzero(np.all(miller1==hkl,axis=-1)))
                print(f"2 num points for {hkl}",np.count_nonzero(np.all(miller2==hkl,axis=-1)))
            I1[i] = np.average(I1_all[np.all(miller1== hkl,axis=-1)])
            if compare_dir!=None:
                I2[i] = np.average(I2_all[np.all(miller2== hkl,axis=-1)])
        all_I_real = I1
        all_I_ideal = I2

        if not get_R_only and results_for_plotting:
            add_screen_properties()
            plt.show()  

            
        if compare_handle != None:

            if not get_R_only and results_for_plotting:
                #print("Plotting orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
                #plot_sectors(sector_histogram h= sector_histogram)
                non_zero_denon_histogram = sector_den_histogram.copy()
                non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
                sector_histogram = np.divide(sector_num_histogram,non_zero_denon_histogram)
                print("Plotting total R factor (incorrect, need to normalise I's)")
                plot_sectors(sector_histogram = sector_histogram)

                #print("plotting full ring orientation-averaged R factor") # Doesn't work because when sector is empty it reduces the average.
                #plot_sectors(R_histogram)

                
                non_zero_denon_histogram = R_den_histogram.copy()
                non_zero_denon_histogram[non_zero_denon_histogram== 0] = 1        
                R_histogram = np.divide(R_num_histogram,non_zero_denon_histogram)       
                print("plotting full ring total R factor (incorrect, need to normalise I's)") 
                plot_sectors(R_histogram)   

            if neutze_R:
                sqrt_ideal = np.sqrt(all_I_ideal)
                sqrt_real = np.sqrt(all_I_real)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                # print("sqrt_ideal",sqrt_ideal)
                # print("-------")
                # print("-------")
                # print("sqrt_real normed",sqrt_real*inv_K) 

                print("---------------")
                print("R: ",R)
                print("---------------")

                # print("R:") (not normalised)
                # R = np.sum(np.abs((sqrt_real - sqrt_ideal)))/np.sum(sqrt_ideal)
                # print(R)
                if not get_R_only:
                    print("sum of (sqrt) real","{:e}".format(np.sum(sqrt_real)),"sum of (sqrt) ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(np.sum(np.abs((sqrt_real - sqrt_ideal)))))
                    #neutze_histogram = np.histogram2d(phi, radial_axis, weights=R, bins=(np.array([-np.pi,np.pi]), radial_edges))[0]             
                    #plot_sectors(neutze_histogram)

                # Pearson Correlation Coefficient:
                x = all_I_ideal
                y = all_I_real
                x_bar = np.mean(x)
                y_bar = np.mean(y)   
                num = np.sum((x-x_bar)*(y-y_bar))
                den = np.sqrt(np.sum(x-x_bar)**2*np.sum(y-y_bar)**2)
                cc = num/den

                return R,cc,None # resolution lims not implemented
    
    ## Continuous (SPI)
    else:
        if SPI_result1 is None:
            print("error, expected SPI but SPI_result1 is None")
            return     
        result1 = copy.deepcopy(SPI_result1)
        result2 = copy.deepcopy(SPI_result2)
        if normalize_to_centre:
            # # Normalize I to 1 at centre 
            # # Average out centre pixels/select central pixels depending on axis dimension's integer parity.
            # A,B,C,D = ( int(np.floor((result1.I.shape[0]+1)/2)), int(np.ceil((result1.I.shape[0]+1)/2)), 
            #             int(np.floor((result1.I.shape[1]+1)/2)), int(np.ceil((result1.I.shape[1]+1)/2)) )
            # result1.I/=np.average(result1.I[A:B+1, C:D+1])
            # if result2 != None:
            #     result2.I/=np.average(result2.I[A:B+1, C:D+1])  
            result1.I/=result1.zero_angle_I
            result1.zero_angle_I = 1
            if result2 is not None:
                result2.I/=result2.zero_angle_I
                result2.zero_angle_I = 1
        else:
            # Normalize  total magnitude of intensity to 1.
            result1.I/=np.sum(result1.I)
            result1.zero_angle_I/=np.sum(result1.I)
            if result2 is not None:
                result2.I/=np.sum(result2.I)
                result2.zero_angle_I/=np.sum(result2.I)
        ## Square grid
        if len(result1.I.shape) > 1:#if type(result1) == Results_SPI:
            if spi_full_rings_only:
                # for calculating R we remove the non-full rings - we represent this visually:
                result1.I *=  (result1.full_ring_mask + 0.1)/1.1        
                if result2 is not None:
                    result2.I *= (result2.full_ring_mask+0.1)/1.1 

            if log_I: 
                z1 = log_function(result1.I)   
                result1.zero_angle_I = log_function(result1.zero_angle_I)
                z1[result1.I == 0] = None
                if result2 != None:
                    z2 = log_function(result2.I)
                    result2.zero_angle_I = log_function(result2.zero_angle_I)
                    z2[result2.I == 0] = None
            #else:
                #TODO fix up this masked stuff i didnt finish 
                # z1 -= cutoff_log_intensity
                # z1[z1<0] = 0
                #z1 = result1.I.copy()/np.sum(result1.I)    
                #if result2 != None:
                    #z2 -= cutoff_log_intensity
                    #z2[z2<0] = 0
                    #z2 = result2.I.copy()/np.sum(result2.I)        
            if log_I and cutoff_log_intensity != None:
                #z1 = np.ma.array(z1, mask=(z1<cutoff_log_intensity)*(np.isnan(z1)))
                z1 = np.ma.array(z1, mask=((z2<cutoff_log_intensity)|(np.isnan(z2))))
                if result2 != None:
                    #z2 = np.ma.array(z2, mask=(z2<cutoff_log_intensity)*(np.isnan(z2)))
                    z2 = np.ma.array(z2, mask=((z2<cutoff_log_intensity)|(np.isnan(z2))))
            else:
                z1 = np.ma.array(z1,mask = np.isnan(z1))
                if result2 is not None:
                    z2 = np.ma.array(z2,mask = np.isnan(z2)) 

            if spi_full_rings_only:
                # Now remove the non-full rings
                result1.I *=  result1.full_ring_mask     
                if result2 is not None:     
                    result2.I *= result2.full_ring_mask
                        
            if result2 != None:
                alpha2 = (result2.full_ring_mask + 1)/2             
            print("Result 1 (Damaged):")
            print("Total screen-incident intensity:","{:e}".format(np.sum(result1.I)))

            if result2 != None:
                combined_data = np.array([z1,z2])
            else: 
                combined_data = z1
            z_min, z_max = np.nanmin(combined_data), np.nanmax(combined_data)      
            #print("DEBUG",z_max,log_range,z_min)
            if log_I:
                z_min = max(z_max-log_range,z_min)
                
            current_cmap = plt.colormaps.get_cmap("plasma")
            current_cmap.set_bad(color='black')
            if not get_R_only:
                #print(f"vmin, vmax: {z_min}, {result1.zero_angle_I}")
                z1_map = plt.imshow(z1,vmin=z_min,vmax=result1.zero_angle_I,cmap=current_cmap)

                if show_grid: # Resolution ring annotations
                    radius = len(result1.q)/2
                    skip = int(radius/5)
                    grid_ring_radii = list(reversed(range(1,int(np.floor(radius))+1)))[::skip]
                    print(grid_ring_radii)
                    for r in grid_ring_radii:
                        if len(result1.q)%2==1:
                            resol = q_to_res(result1.q[int(radius)+r][int(radius)])*ang_per_bohr
                            midpoint=radius-1
                        else: # rough
                            midpoint=radius-0.5
                            resol = 0
                            for x in (int(np.floor(radius)),int(np.ceil(radius))):
                                row=x+r
                                if x+r >= len(result1.q):
                                    row-=1
                                resol += 0.5 * q_to_res(result1.q[row][x])*ang_per_bohr
                        theta=np.linspace(0, 2*np.pi, 1000)
                        label=f"${resol:.1f} Å$"
                        plt.gca().plot(midpoint + r*np.cos(theta), midpoint + r*np.sin(theta),
                                       alpha = 0.8,
                                       dashes=(5,10),
                                       color="silver",label=label) #whitesmoke lightgray silver
                        plt.gca().text(midpoint,midpoint+r,label,color="w",alpha=1,
                                       horizontalalignment='center',
                                       verticalalignment='top')
                        plt.gca().set_ylim(-0.5,result1.q.shape[0]-0.5)
                        plt.gca().set_xlim(-0.5,result1.q.shape[1]-0.5)
                    #labellines.labelLines(plt.gca().get_lines())
                
                plt.colorbar(z1_map, label="Intensity (arb. units)",
                             format=lambda x, _: "$10^{"+f"{x:.0f}"+"}$")
                            
                plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                if plot_handle!="":
                    plt.savefig("I_"+plot_handle+".png",format="png",dpi=dpi)
                    plt.savefig("I_"+plot_handle+".pdf",format="pdf",dpi=dpi)
                else:
                    plt.savefig("tmp_I_real.pdf",format="pdf",dpi=dpi)
                if show_plot:
                    plt.show()
            if result2 != None:       
                if not get_R_only:
                    print("Result 2 (Undamaged):")
                    print("Total screen-incident intensity:","{:e}".format(np.sum(result2.I)))
                    z2_map = plt.imshow(z2,vmin=z_min,vmax=result2.zero_angle_I,cmap=current_cmap)
                    plt.colorbar(z2_map)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_I_ideal.pdf",format="pdf",dpi=dpi)
                    plt.show()
                    print("R:")
                    fig, ax = plt.subplots()
                    I_tmp = result2.I
                    if log_I:
                        I_tmp = z2.copy()
                    alpha = (I_tmp-np.nanmin(I_tmp))/(np.nanmax(I_tmp) - np.nanmin(I_tmp))
                    alpha[np.isnan(alpha)] = 0
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # Scales I_real to I_ideal's tot intensity
                R_cells = np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal))
                R = np.sum(R_cells)
                print(R)               
                if not get_R_only:
                    #min_R_dmg_pixel = 0.1 # The minimum R factor contributed by a pixel (multiplied by the number of pixels), for it to be displayed.
                    alpha_prop_to_I = False
                    bg = np.full((*z1.shape, 3), 0, dtype=np.uint8) #bg = np.full((*z1.shape, 3), 70, dtype=np.uint8)

                    ax.imshow(bg)
                    R_cells *= len(R_cells)**2# multiply by num cells to give the 'weighted contribution' (such that R is now like the weighted average) 
                    alpha[R_cells < min_R_dmg_pixel] = 0  # hide insignificant pixels
                    if not alpha_prop_to_I:
                        alpha[alpha != 0] = 1  # override.
                    R_map = ax.imshow(R_cells,vmin=0,vmax=0.4,alpha=alpha,cmap=cmap)     
                    plt.colorbar(R_map)
                    ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                    ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                    ticks = ticks[len(ticks)//10::len(ticks)//5]
                    ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]                    
                    plt.xticks(ticks,ticklabels)
                    plt.yticks(ticks,ticklabels)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_R_contributions.pdf",format="pdf",dpi=dpi)
                    plt.show()

                    R_num = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)))        

                    norm_root_diff_map = True
                    if norm_root_diff_map:  #TODO clean mess
                        print("Plotting normalised root difference map")
                        fig, ax = plt.subplots()
                        R_cells = np.abs((inv_K*sqrt_real - sqrt_ideal)/sqrt_ideal)
                        ax.imshow(bg)
                        R_map = ax.imshow(R_cells,vmin=0,vmax=2,alpha=alpha,cmap=cmap2)     
                        plt.colorbar(R_map)
                        ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                        ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                        ticks = ticks[len(ticks)//10::len(ticks)//5]
                        ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]                        
                        plt.xticks(ticks,ticklabels)
                        plt.yticks(ticks,ticklabels)
                        plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                        plt.savefig("tmp_R_pixel.pdf",format="pdf",dpi=dpi)
                        plt.show()   
                        # print("Plotting change map") 
                        # fig, ax = plt.subplots()
                        # R_cells = inv_K*sqrt_real/sqrt_ideal
                        # ax.imshow(bg)
                        # R_map = ax.imshow(R_cells,vmin=None,vmax=2,alpha=alpha,cmap=cmap2)     
                        # plt.colorbar(R_map)
                        # ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                        # ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                        # plt.xticks(ticks,ticklabels)
                        # plt.yticks(ticks,ticklabels)
                        # plt.show()                          
                    print ("Plotting log ratio")
                    I1 = result1.I #real
                    I2 = result2.I #ideal
                    # Average out centre pixels/select central pixels depending on axis dimension's integer parity.
                    A,B,C,D = ( int(np.floor((I1.shape[0]+1)/2)), int(np.ceil((I1.shape[0]+1)/2)), 
                                int(np.floor((I1.shape[1]+1)/2)), int(np.ceil((I1.shape[1]+1)/2)) )
                    I1_divisor = np.average(I1[A:B+1, C:D+1])
                    I2_divisor = np.average(I2[A:B+1, C:D+1])                    
                        
                    log_ratio = log_function((I1/I1_divisor)/(I2/I2_divisor))
                    fig, ax = plt.subplots()
                    ax.imshow(bg)
                    norm = TwoSlopeNorm(vmin=log_diff_vmin, vcenter=0, vmax=log_diff_vmax)
                    I_map = ax.imshow(log_ratio,norm=norm,cmap = cmap2)#cmap="plasma")#"nipy_spectral_r")
                    cb = plt.colorbar(I_map)
                    cb.ax.set_yscale('linear')
                    ticks = np.linspace(0,len(result1.xy)-1,len(result1.xy))
                    ticklabels = ["{:6.2f}".format(q_row_0_el[1]) for q_row_0_el in result1.q_xy[0]]
                    ticks = ticks[len(ticks)//10::len(ticks)//5]
                    ticklabels = ticklabels[len(ticklabels)//10::len(ticklabels)//5]
                    plt.xticks(ticks,ticklabels)
                    plt.yticks(ticks,ticklabels)
                    plt.gcf().set_figwidth(custom_fig_width); plt.gcf().set_figheight(custom_fig_height)
                    plt.savefig("tmp_log_diff.pdf",format="pdf",dpi=dpi)
                    plt.show()                      


                    print("sum of real","{:e}".format(np.sum(sqrt_real)),"sum of ideal","{:e}".format(np.sum(sqrt_ideal)),"sum of abs difference","{:e}".format(R_num))       
                

                ### Calculate damage measures
                # Iterate through bins of resolutions
                resolutions = q_to_res(np.sqrt(np.apply_along_axis(np.sum,2, result1.q_xy**2))) 
                #min_res = np.min(resolutions)  # THIS IS NOT FROM THE RIM!!! It includes all.
                min_res = np.max([np.max(resolutions[0,:]),np.max(resolutions[:,0])])  # rim resolution
                max_res = np.min([6,np.max(resolutions)])
                res_lims = np.linspace(min_res,max_res,20)
                cc = np.zeros(len(res_lims)) # Pearson cc
                R_vect = cc.copy() # R_dmg (for full dataset up to the given resolution)
                binned_R_vect = cc.copy() # R_dmg for individual resolution bins/rings
                binned_q_R_dict = {}
                binned_q_I_ratio_dict  = {}
                # Generate dictionary of evenly spaced q, spanning the range of 0 to q corresponding to min_res. Bins centred on 1/2*delta_q, 3/2*delta_q and so on.
                delta_q = 0
                for i,val in enumerate(np.linspace(0,res_to_q(min_res),20,endpoint=False)):
                    binned_q_R_dict[val+0.5*delta_q] = np.nan
                    if i == 1:
                        delta_q = val
                        binned_q_R_dict.pop(0)
                        binned_q_R_dict.pop(val)
                        binned_q_R_dict[0.5*delta_q],binned_q_R_dict[1.5*delta_q] = [np.nan]*2
                binned_q_I_ratio_dict  = copy.deepcopy(binned_q_R_dict)                        
                
                for i in range(len(res_lims)):
                    lower = resolutions * 0 + res_lims[i]
                    upper = np.Infinity
                    # Select data up to some resolution limit
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    # cc
                    x_bar = np.mean(x)
                    y_bar = np.mean(y)   
                    num = np.sum((x-x_bar)*(y-y_bar))
                    den = np.sqrt(np.sum((x-x_bar)**2)*np.sum((y-y_bar)**2))
                    cc[i] = None
                    if den != 0:
                        cc[i] = num/den
                    # R factor
                    sqrt_real = np.sqrt(x)
                    sqrt_ideal = np.sqrt(y)
                    inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    R_vect[i] = None
                    if np.sum(sqrt_ideal) != 0:
                        R_vect[i] = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))                    
                    # R factor but binned to resolution (note the resolutions still correspond to the maximum scattering angles for the bin, i.e. the bin is not centred on the corresponding q for the resolution.)
                    delta_res = res_lims[1]-res_lims[0]
                    lower = resolutions * 0 + res_lims[i]
                    upper = resolutions * 0 + res_lims[i] + delta_res/2
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    binned_sqrt_real = np.sqrt(x)
                    binned_sqrt_ideal = np.sqrt(y)
                    binned_inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    binned_R_vect[i] = None                    
                    if np.sum(binned_sqrt_ideal) != 0:
                        binned_R_vect[i] = np.sum(np.abs((binned_inv_K*binned_sqrt_real - binned_sqrt_ideal)/np.sum(binned_sqrt_ideal)))
                # R factor but bins separated by constant delta q.
                for key in binned_q_R_dict.keys():
                    lower = resolutions*0 + q_to_res(key+delta_q/2)
                    upper = resolutions*0 + q_to_res(key-delta_q/2)
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    #tmp = copy.deepcopy(result1.I)
                    #tmp[tmp == 0] = np.nan
                    #plt.imshow(tmp)
                    sqrt_real = np.sqrt(x)
                    sqrt_ideal = np.sqrt(y)
                    inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real)   # normalises I_real to I_ideal's tot intensity
                    binned_q_R_dict[key] = None   
                    if np.sum(sqrt_ideal) != 0:
                        binned_q_R_dict[key] = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))                         
                # Average intensity
                I1_centre = result1.I[int(len(result1.I)/2)][int(len(result1.I)/2)]
                I2_centre = result2.I[int(len(result2.I)/2)][int(len(result2.I)/2)]
                for key in binned_q_I_ratio_dict.keys():
                    lower = resolutions*0 + q_to_res(key+delta_q/2)
                    upper = resolutions*0 + q_to_res(key-delta_q/2)
                    x = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result1.I)
                    y = np.extract((resolutions >= lower) * (resolutions < upper)*result1.I*result2.I >np.zeros(result1.I.shape),result2.I)
                    if np.sum(y) != 0:
                        binned_q_I_ratio_dict[key] = log_function(np.mean( (x/I1_centre)/(y/I2_centre) ))                    

                #plt.plot((bin_lims[0:-1] + bin_lims[1:])/2,cc)
                if not get_R_only:
                    plt.rcParams["text.usetex"] = True
                    plt.plot(res_lims,R_vect,color="red")
                    plt.ylabel("R")
                    plt.xlabel("d ($\AA$)")                    
                    plt.ylim(0,1.04)
                    plt.gca().invert_xaxis()
                    plt.show()                    
                    plt.plot(res_lims,cc)
                    plt.ylabel("CC")
                    plt.xlabel("d ($\AA$)")
                    plt.ylim(0,0.2)
                    plt.gca().invert_xaxis()
                    plt.show()
                damage_dict  = dict(
                    R = R_vect,
                    bin_res_R = binned_R_vect,
                    cc = cc,
                    resolutions = res_lims,
                    bin_q_R = binned_q_R_dict,
                    bin_q_I = binned_q_I_ratio_dict
                )
                return damage_dict
            
        ## Circle grid
        else:
            if log_I: 
                z = log_function(result1.I)
            else:
                z = result1.I        
            if log_I and cutoff_log_intensity != None:
                #cutoff_log_intensity = -1
                z -= cutoff_log_intensity
                z[z<0] = 0
                pass
            radial_axis = result1.X
            if plot_against_q:
                radial_axis = result1.q         
            fig = plt.figure()
            ax = fig.add_subplot(projection="polar")
            # phi_mesh = result.phi_mesh
            # if crystal_aligned_frame:
            #     phi_mesh = result.phi_aligned_mesh   
            #     print("crystal aligned not implemented for SPI yet")    
            ax.pcolormesh(result1.phi, radial_axis, z,cmap=cmap)
            ax.plot(result1.phi, radial_axis, color = 'k',ls='none')
            plt.grid(alpha=min(show_grid,0.6),dashes=(5,10)) 

            if radial_lim:
                bottom,top = plt.ylim()
                plt.ylim(bottom,radial_lim)
            if log_radial:
                plt.yscale("log")

            if result2 != None:
                print("R:")
                sqrt_real = np.sqrt(result1.I)
                sqrt_ideal = np.sqrt(result2.I)
                inv_K = np.sum(sqrt_ideal)/np.sum(sqrt_real) 
                R = np.sum(np.abs((inv_K*sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                print(R)   

                # Pearson Correlation Coefficient:
                x = result1.I
                y = result2.I
                x_bar = np.mean(x)
                y_bar = np.mean(y)   
                num = np.sum((x-x_bar)*(y-y_bar))
                den = np.sqrt(np.sum(x-x_bar)**2*np.sum(y-y_bar)**2)
                cc = num/den                             
                return R,cc,None  # resolution lims not implemented
                # print("R: (not normalised)")
                # R = np.sum(np.abs((sqrt_real - sqrt_ideal)/np.sum(sqrt_ideal)))
                # print(R)

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
    matplotlib.colormaps.register(cmap=newcmap)

    return newcmap

# def get_result(filename,results_dir,compare_dir = None):
#     #Requires all orientations of result_handle in compare_handle, but not vice versa.
#     fpath = os.path.join(results_dir, filename)
#     if os.path.isfile(fpath):
#         with open(fpath,'rb') as f:
#             result1 = pickle.load(f)
#     else: 
#         return "__PASS__" "__PASS__"
#     result2 = None
#     if compare_dir != None:
#         if filename in os.listdir(compare_dir):
#             fpath2 = os.path.join(compare_dir, filename)
#             with open(fpath2,'rb') as f:
#                 result2 = pickle.load(f)     
#                 result1.diff(result2)             
#         else:
#             print("ERROR, missing matching orientation in comparison directory")
#             return None, None # No corresponding file found.          
#     return result1, result2

def get_result(filename,results_dir,compare_dir = None):
    return Results.get_result(filename,results_dir,compare_dir)


##### https://scripts.iucr.org/cgi-bin/paper?S0021889807029238, http://superflip.fzu.cz/

import pandas as pd
import csv
def create_reflection_file(result_handle,results_parent_dir = RESULTS_LOCAL_PATH,out_directory="reflections/",overwrite=False,artificial_I_scale=1,symmetry_override=None):
    '''
    Generates a .rfl file, compatible with Superflip
    '''
    print("Creating reflection file for",result_handle)
    results_dir =  path.abspath(path.join(__file__ ,"../")) + "/"+ results_parent_dir + result_handle+"/"
    assert path.isdir(results_dir), f"Directory not found: {results_dir}" 
    out_directory = path.abspath(path.join(__file__ ,"../")) + "/"+ out_directory
    os.makedirs(out_directory, exist_ok=True) 
    init = False
    for filename in os.listdir(results_dir):
        result:Results = get_result(filename,results_dir)[0]
        if result == "__PASS__":
            continue
        if result == None:
            break          
        miller_indices = result.miller_indices
        intensity = np.array([result.I]).T * artificial_I_scale
        print("Shape I",intensity.shape)
        data = np.concatenate((miller_indices,intensity),axis=1) #TODO work with separate times
        # Put in data frame
        columns = ["h","k","l","I"]
        new_df = pd.DataFrame(data=data, columns=columns)
        if init == False:
            df = new_df
            init = True
        else:
            df = pd.concat((df,new_df))
    # Save as file
    out_path = out_directory + result_handle
    if os.path.isfile(out_path): 
        if not overwrite:
            print("Cannot write, file already present at",out_path)
            return
        os.remove(out_path)
    #columns.reverse() #???
    df = df.sort_values(by=["l","k","h"],axis=0)
    for i in ("hkl"):
        df[i] = df[i].astype('int')
    df["I"] = df["I"].astype('float')
    df = df.round(6)
    #df.drop_duplicates(subset = ["h","k","l"],inplace=True) # TODO should take average.
    #df = df[df['I']>=0.01] # Update: This is silly TODO temporary fix for appearance of low values that needs to be squashed.
    
    cell_geom = []
    symmetry = result.symmetry if symmetry_override is None else symmetry_override
    for l in (result.cell_dims, result.cell_angles,symmetry):
        cell_geom += [str(v) for v in l]
    cell_geom = ' '.join(cell_geom)
    print(cell_geom)

    def make_file(_path,_df:pd.DataFrame):
        with open(_path,'w') as f:
            f.write(f'# {cell_geom}\n')
        _df.to_csv(_path,mode='a',header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")
        
    make_file(out_path+"_unmerged.rfl",df)
    #df.to_csv(,mode='a',header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")
    

    df_merged = df.groupby(["h","k","l"]).mean().reset_index()
    print("shape merged",df_merged.shape)
    df_merged = df_merged.sort_values(by=["l","k","h"],axis=0)
    #df_merged.to_csv(out_path+".rfl",mode='a',header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")
    make_file(out_path+".rfl",df_merged)

class ScalingMethod:
    def __init__():
        pass
    def get_I_norm(self,rfl_file_path,photon_energy="same"):
        print("implement")
        raise Exception()
    def get_scaled(self,h,k,l):
        print("implement")
        raise Exception()


def noise_and_photon_count_curve(x,a,b,c,d):
    if a < 0:  # Enforce positive photon counting error. Probably not how you're meant to do it.
        return -9999999
    if b < 0 or d < 0:  # Enforce that other sources of error are positive
        return -9999999
    y = a * x**0.5 +b*x**c + d 
    return y


class ScalingByReference(ScalingMethod):
    '''
    Scales I based on reference reflection photon count.
    '''
    # def __init__(self,h,k,l,I_meas,I_sigma,photon_energy):
    #     self.reference_reflection = (h,k,l,I_meas,I_sigma)
    #     self.photon_energy = photon_energy
        
    #     self.photon_count = 
    def __init__(self,h,k,l,I_meas,I_sigma,photon_energy):
        self.reference_reflection = (h,k,l,I_meas,I_sigma)
        self.photon_energy = photon_energy
        
        pass

    def get_I_norm(self,rfl_file_path,photon_energy="same"):
        '''
        A .rfl file only ever has one set of reflections, either from single experiement or merged data. Regardless, we make the same treatment. 
        '''
        if photon_energy == "same":
            photon_energy = self.photon_energy
        
        #photon_count = I/

def read_rfl_file(rflFile):
    data_dict = dict (h=0,
                    k=1,
                    l=2,
                    I_sim=3,
                )
    rows =[]
    with open(rflFile, 'r') as f:
        for line in f:
            data = line.split()
            row = []
            for k,idx in data_dict.items():
                row.append(float(data[idx]))
            rows.append(row)
    return pd.DataFrame(rows, columns=list(data_dict.keys()))

class ScalingByCopyingSigmaRatio(ScalingMethod): # Should be valid for ideal sim. 
    def __init__(self,cifFile,rflFile):
        df_sim = read_rfl_file(rflFile)#.astype(int)
        df_real = read_cif(cifFile)#.astype(int)
        #self.df = pd.concat([df_sim, df_real], ignore_index=True, sort=False)
        self.df = pd.merge(df_sim,df_real, on=['h','k','l'])
        assert self.df.shape[1]==6, f"{self.df.shape}"

    def get_scaled(self,h,k,l):
        df = self.df

        data = df[df['h']==float(h)][df['k']==float(k)][df['l']==float(l)]
        # if not (data['I_meas'].notnull().all() and data['I_sigma'].notnull().all()):
        #     return None,None
        if data.shape[0]==0:
            return None,None
        I_scaled = data['I_sim']*self.get_I_norm()
        I_sigma = data['I_sigma']/data['I_meas']*I_scaled
        assert I_scaled.shape[0] == I_sigma.shape[0] == 1
        I_scaled = I_scaled.iloc[0]
        I_sigma= I_sigma.iloc[0]
        assert I_scaled==I_scaled
        assert I_sigma==I_sigma
        return I_scaled,I_sigma
    def get_I_norm(self):
        # compare total irradiance over shared points.
        df = self.df[self.df['I_meas'].notnull()][self.df['I_sim'].notnull()]
        return df['I_meas'].sum()/df['I_sim'].sum()



# Won't work... we need different fluences
# class ScalingByFit(ScalingMethod):
#     def __init__(self):
#         self.curve=noise_and_photon_count_curve
#         self.popt=None
#         self.pcov=None
#     def get_I_norm(self,rfl_file_path,photon_energy="same"):
#         self.curve_fit(cifFile="TODO")

#         total_I, reflections_observed = getTotalReflectionI(cifFile)
#         I_norm = total_I/
#         return number 
#     def curve_fit(self,cifFile="TODO"):
#         #TODO IMPLEMENT
#         self.popt,self.pcov = (1.42845205e-06, 2.78574265e-01, 9.21457720e-01, 1.91285876e+01),None
#         self.plot_curve_fit()
#         #scipy.optimize.curve_fit(curve,df['I_meas'],df['I_sigma'],sigma=df['I_meas'], absolute_sigma=True)
#     def plot_curve_fit(self,log=False):
#         if log:
#             x = np.logspace(0,5,100)
#         else:
#             x = np.linspace(0,5e4,100)
#         plt.scatter(x,self.curve(x,*self.popt),s=1)


def rfl_to_sca(result_handle, reflections_dir = "reflections/", out_directory = "scalepack/",overwrite=True,detector_gain=0.6,scaling_method : ScalingMethod = None,create_mtz=True):
    '''Converts .rfl file to scalepack .sca file
    See https://www.ccp4.ac.uk/html/scala.html#files
    '''
    out_directory=path.abspath(path.join(__file__ ,"../")) + "/"+ out_directory #XXX change default to none. assume not none is an absolute path
    reflections_dir=path.abspath(path.join(__file__ ,"../")) + "/"+ reflections_dir #XXX
    os.makedirs(out_directory, exist_ok=True) 
    rfl_file_path = reflections_dir + result_handle + ".rfl"
    assert path.isfile(rfl_file_path), f"Reflection file not found: {rfl_file_path}" 

    out_path = out_directory + result_handle + ".sca"
    save_action = "x"
    if overwrite:
        save_action = "w"    
    # First find max length of intensities, so can scale values down.
    max_length = 0
    with open(rfl_file_path, 'r') as f_a:
        for line in f_a:
            entries = line.split()
            assert len(entries) >=3, entries
            if entries[0]+entries[1]+entries[2] == '0'*3: 
                continue # Ignore (0,0,0) reflection
            
            max_length = max(max_length,len(line.split()[3].split('.')[0]))
    #
    #I_norm = scaling_method.get_I_norm(rfl_file_path)
    # with open(rfl_file_path, 'r') as f_a:
    #     for i, line in enumerate(f_a):
    #         print(i,line)
    with open(rfl_file_path, 'r') as f_a, open(out_path, save_action) as f_b:
        # placeholder boilerplate 
        indent = ' '*3
        f_b.write(indent+' 1\n -987\n')
                
        #CRYST1   79.200   79.200   37.800  90.00  90.00  90.00 P 1           1
        
        #b.write(indent+' 79.000    79.000    38.000    90.000    90.000    90.000 P43212\n')
        #b.write(indent+' 79.000    79.000    38.000    90.000    90.000    90.000 P1\n')
        # Miller indices (hkl) | IMEAN_dataset | SIGIMEAN_dataset
        for i, line in enumerate(f_a):
            if i == 0:
                entries = line.split()[1:]  # format '# a b c A B C P x [x] [x]'
                a,b,c,A,B,C = [float(s) for s in entries[:6]]
                symm = ''.join(entries[6:])
                # Cell geometry TODO integrate with actual input.]
                def F(_c):
                    return str(f"{_c:.6f}")[:6]
                f_b.write(indent + f' {F(a)}    {F(b)}    {F(c)}    {F(A)}    {F(B)}    {F(C)} {symm}\n')
                continue

            # Initialise elements
            h=k=l = ' '*4
            I_mean=sigI_mean = ' '*8
            

            # convert data to usable format from .rfl file
            entries = line.split()
            #   hkl
            if entries[0]+entries[1]+entries[2] == '0'*3:
                continue # Ignore (0,0,0) reflection
            #   I
            if scaling_method is not None:
                Isim, Isigma = scaling_method.get_scaled(*entries[:3])
                if Isim is None:
                    continue
            else:
                I_scaling_power = max_length - 7
                Isim = float(entries[3])/10**I_scaling_power
                Isigma = np.sqrt(Isim)

            entries[3] ='%.0f'%(Isim)
            # I sigma
            entries.append('%.1f'%Isigma) # .rfl doesn't have sigma.

            


            # Populate elements
            for i,q in enumerate([h,k,l,I_mean,sigI_mean]):
                assert len(entries[i]) < len(q), "Error, entry of '"+entries[i]+"' has length >= "+str(len(q))   # Use '<' not '<=' because need a space.
                q = q[:-len(entries[i])]+ entries[i]
                f_b.write(q)
            f_b.write('\n')
    
    
    # Create mtz file too using phenix (if installed).
    mtz_file = None
    if create_mtz:
        mtz_file = f"{path.abspath(out_path)[:-4]}.mtz"
        create_mtz_args =[
            "phenix.reflection_file_converter",
            path.abspath(out_path),
            f"--mtz={mtz_file}",
        ]
        print (f"Running: {' '.join(create_mtz_args)}")
        subprocess.run(create_mtz_args,stdout=subprocess.PIPE)

    return path.abspath(out_path), mtz_file

# For some reason passing in the miller indices does not produce the right map, probably due to wrong phases since R factor is unchanged.
# So need to specify resolution sadly. 
# Unfortunately this also makes working with it slower.
#def phenix_fcalc(pdb_file,miller_indices_mtz):
def phenix_fcalc(pdb_file,high_resolution,real=False):
    type_tag = "xyz" if real else "cplx" 
    out_file_name = f"{pdb_file[:-4]}_fcalc_{type_tag}.mtz"
    if path.exists(out_file_name):
        os.remove(out_file_name)
    args =[
        "phenix.fmodel",
        pdb_file,
        #miller_indices_mtz,
        f"high_resolution={high_resolution}",
        f"file_name={out_file_name}",
        "use_asu_masks=False",
        "algorithm=direct",
        "type=real" if real else "type=complex",
        "obs_type=amplitudes",
        "grid_resolution_factor=1/3",
        "ignore_hydrogens=False",
        #"data_column_label=FOBS,SIGFOBS",
        #f"high_res={high_resolution}",
    ]
    print (f"Running: {' '.join(args)}")
    subprocess.run(args,stdout=subprocess.PIPE)#,stdout=log)
    return out_file_name

def phenix_fcalc_from_file(pdb_file,miller_indices_mtz,real=False):
    assert real
    type_tag = "xyz" if real else "cplx" 
    out_file_name = f"{pdb_file[:-4]}_fcalc_{type_tag}.mtz"
    if path.exists(out_file_name):
        os.remove(out_file_name)
    args =[
        "phenix.fmodel",
        pdb_file,
        miller_indices_mtz,
        f"file_name={out_file_name}",
        "use_asu_masks=False",
        "algorithm=direct", # *fft direct
        "type=real" if real else "type=complex",
        "obs_type=amplitudes",
        "grid_resolution_factor=1/3",
        "ignore_hydrogens=False",
        #"data_column_label=FOBS,SIGFOBS",
    ]
    print (f"Running: {' '.join(args)}")
    subprocess.run(args,stdout=subprocess.PIPE)
    return out_file_name



def phenix_R(pdb_file,reflections):
    args =[
        "phenix.model_vs_data",
        pdb_file,
        reflections,
        #f"high_res={high_resolution}",
    ]
    print (f"Running: {' '.join(args)}")
    proc = subprocess.run(args,encoding='utf-8',stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    found_line = False
    for line in proc.stdout.split('\n'):
        if line.startswith("  r_work:"):
            print(line)
            found_line=True
    if not found_line:
        #print(proc.stdout)
        print("Error! Are occupancies zero? Are symmetries the same?")

def read_scalepack(result_handle,scalepack_dir = "scalepack/",skip_header=3):
    file_path = path.abspath(path.join(__file__ ,"../")) + "/"+ scalepack_dir + result_handle + ".sca"
    indices = []
    I = []
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if i < skip_header:
                continue
            indices.append([
                int(line[0:4]),
                int(line[4:8]),
                int(line[8:12])
            ])
            I.append(float(line[12:20]))
            
    return np.array(indices), np.array(I)

def read_hkl(fpath):
    ''' Read h k l from shelxl hkl-format file '''
    indices = []
    with open(fpath, 'r') as f:
        for line in f:
            h = int(line[0:4])
            k = int(line[4:8])            
            l = int(line[8:12])
            indices.append([h,k,l])            
    return np.array(indices)

#create_reflection_file("hen_v7__eal",True)
#rfl_to_sca("hen_v7_real")

#####
# stylin' 
def stylin(exp_name1,exp_name2,radial_lim,damaged_and_undamaged=False,get_R_only = False,SPI=False,SPI_max_q=None,SPI_result1=None,SPI_result2=None,results_parent_dir = RESULTS_LOCAL_PATH,show_labels=False,cutoff_log_intensity = None,**kwargs):
    experiment1_name = exp_name1#"Lys_9.95_random"#exp_name1
    experiment2_name = exp_name2#"lys_9.80_random"#exp_name2 

    #####
    results_dir = path.abspath(path.join(__file__ ,"../")) + "/"+ results_parent_dir +experiment1_name+"/"
    for filename in os.listdir(results_dir):
        result1,result2 = get_result(filename,results_dir)
        if not (result1.for_plotting):
            get_R_only = True
            break

    font = {'family': 'serif',
            'size'   : 10}

    plt.rc('font', **font)

    use_q = False # TODO remove this option
    log_radial = False
    log_I = True
    #cutoff_log_intensity = -1#-1
    
    cmap = ""
    if not get_R_only:
        try:
            cmap = shiftedColorMap(matplotlib.cm.RdYlGn_r,midpoint=0.2,name="shiftedcmap")#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
        except Exception as e:
            try:
                cmap =  plt.get_cmap("shiftedcmap")
            except Exception as e2:
                print(e)
                print(e2)
                cmap =  plt.get_cmap("RdYlGn_r")
        cmap.set_bad(color='black')
    cmap_power = 1.6
    min_alpha = 0.3
    max_alpha = 1
    colour = "y"
    full_crange_sectors = False

    cmap_intensity = "inferno"


    # screen_radius = 150#55#165    #
    # q_scr_lim = experiment.r_to_q(screen_radius) #experiment.r_to_q_scr(screen_radius)#3.9  #NOTE this won't be the actual max q_parr_screen but ah well.
    # zoom_to_fit = True
    # ####### n'
    # # plottin'

    # as q = ksin(theta).  
    #experiment1.max_q*(1/np.sqrt(2)) # for flat screen. as q_scr = qcos(theta). (q_z = qsin(theta) =ksin^2(theta)), max theta is 45. (Though experimentally ~ 22 as of HR paper)

    zoom_to_fit = False
    if not zoom_to_fit:
        radial_lim*=1.02
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
    exp_1_tag = f" ({exp_name1})"
    exp_2_tag = f" ({exp_name2})"
    if not SPI:
        if not get_R_only:
            print("----R Sectors unaligned----")
            scatter_scatter_plot(crystal_aligned_frame = False,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True,results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs) 
            print(f"----Intensity of experiment 1{exp_1_tag}----")
            scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print(f"----Intensity of experiment 2{exp_2_tag}----")
            scatter_scatter_plot(crystal_aligned_frame = False,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment2_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print(f"----Intensity of experiment 1{exp_1_tag} aligned----") 
            scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=False,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print(f"----Intensity of experiment 2{exp_2_tag} aligned----")
            scatter_scatter_plot(crystal_aligned_frame = True,show_grid = True, num_arcs = 25, num_subdivisions = 40,result_handle = experiment2_name, fixed_dot_size = False, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=0.5,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap_intensity,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
            print("----R Sectors aligned----")
            damage_dict = scatter_scatter_plot(crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)
        else:
           print("GET R ONLY",get_R_only)
           damage_dict = scatter_scatter_plot(get_R_only = True,crystal_aligned_frame = True,full_range = full_crange_sectors,num_arcs = 25, num_subdivisions = 40,result_handle = experiment1_name, compare_handle = experiment2_name, fixed_dot_size = True, results_parent_dir=results_parent_dir, cmap_power = cmap_power, min_alpha=min_alpha, max_alpha = max_alpha, solid_colour = colour, crystal_pattern_only = False,show_labels=show_labels,log_dot=True,dot_size=1,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,log_I=log_I,cutoff_log_intensity=cutoff_log_intensity,**kwargs)

    else:
        use_q = False
        log_radial = False
        log_I = True
        #cutoff_log_intensity = None # -1 or None
        cmap = "PiYG_r"# "plasma"
        cmap2 = "seismic"
        radial_lim = SPI_max_q
        damage_dict = scatter_scatter_plot(get_R_only=get_R_only,log_I = log_I, cutoff_log_intensity = cutoff_log_intensity, SPI_result1=SPI_result1,SPI_result2=SPI_result2,radial_lim=radial_lim,plot_against_q = use_q,log_radial=log_radial,cmap=cmap,cmap2=cmap2,**kwargs)


    return damage_dict

def res_to_q(d):
    '''
    Pass d = None to use default resolution (when put into exp_args)
    '''
    if d == None:
        return None
    return 2*np.pi/d
def q_to_res(q):
    return 2*np.pi/q

#TODO 
# log doesnt work atm.
# Get the rings to correspond to actual rings

#%% vvvvv Scattering Playground vvvv
#Scattering Playground

# 0.18081540711890387 SALT
# 0.17781329569346938 NO SALT

if __name__ == "__main__":
    SEEDED = False
    #RANDOM_WATER = True; NUM_RANDOM_WATER = 702;# RANDOM_WATER_EACH_TIME_STEP=True  #TODO implement rand water each time step


    fig_width = 3.49751 # 20
    fig_height = fig_width*3/4 # 20
    ### Simulate
    target_options = ["lys_salt","lys_no_salt","lys_salt_HF","neutze","hen","tetra","glycine","fcc","galliHigh"
                      "lys_nass_probe_35","copper_sulfate"]
    #============------------User params---------==========#

    #R:  0.03453990841341609
    #R:  0.039555455273223315
    #target = "copper_sulfate" #lys_nass_probe_35"#"glycine"  #target_options[2]
    target = "lys_nass_probe_35" #lys_nass_probe_35"#"glycine"  #target_options[2]
    best_resolution = 0.5 # 1.58 (abdullah) # 1.3 # 2   # resolution (determining max q)
    worst_resolution = 30 #None #30 # 'resolution' corresponding to min q

    #### Individual experiment arguments 
    tag = "probe" # Non-SPI i.e. Crystal only, tag to add to folder name. Reflections saved in directory named version_number + target + tag named according to orientation .
    start_time = 0.9 #-18#-18#-12#-6
    end_time = 1.1 #18#12#6
    laser_firing_qwargs = dict(
        # pixel sampling method (Neutze) if True - Miller indices if False
        SPI = False,  # sampling method, if False, bragg spots. if True, detector pixels. TODO change name
        SPI_resolution = best_resolution,
        pixels_across = 100,  # for SPI TODO shld go on xfel exp params.
        # miller
        random_orientation = False, #bragg spot sampling only, TODO refactor to be in same place as other orients...# orientation is synced with second 
    )
    ##### Crystal params
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 "unit" cells per supercell # Bragg spots will be sampled based on the cell scale, not the supercell scale.
        num_supercells = 29*33*37,#100, # 35409
        supercell_simulations = 2, #150
        positional_stdv = 0.1,#0.2,  #Intro   duces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        #supercell_simulations = 900, #150
        #positional_stdv = 0.05,#0.2,  #Intro   duces disorder to positions. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if gauging serial crystallography R factor, as should average out. 0.2 neutze.
        include_symmetries = True,  # should unit cell contain symmetries?
        cell_packing = "SC",
        #rocking_angle = 0.1,  #  (approximating mosaicity - use 0.02 for proper, use a high value, like 1-10, and set a low max triple miller indice to disallow seemingly impossible indices (due to rocking angle/our implementation of it via momentum conservation formulae) that mimic studies that use the first few miller indices )
        rocking_angle = 0.02,  #  (approximating mosaicity - use 0.02 for proper, use a high value, like 1-10, and set a low max triple miller indice to disallow seemingly impossible indices (due to rocking angle/our implementation of it via momentum conservation formulae) that mimic studies that use the first few miller indices )
        #CNO_to_N = True,   # whether the plasma simulation approximated CNO as N  #TODO move this to indiv exp. args or make automatic
        random_waters=None,
        zero_bfactors=False,
        allow_skip_species=True
    )
    crystal_2_has_deviations=False

    show_crystal = True

    #### XFEL params
    #TODO make it so reflections don't overwrite same orientation, as stochastic now.
    energy = 9200#7100 # eV
    exp_qwargs = dict(
        detector_distance_mm = 100,
        screen_type = "flat",#"hemisphere"
        q_minimum = res_to_q(worst_resolution),#None #angstrom
        q_cutoff = res_to_q(best_resolution), #(best_resolution),#2*np.pi/2
        t_fineness=25,   
        #####crystal stuff (miller)
        max_miller_idx = 25, #None, # = m, [overrides max q so given by q with miller indices (m,m,m)]
        all_miller_indices = True, # False, # whether to find all bragg points at or below the max miller index AND between min and max q
        miller_indices_override=None,
        spot_fraction_per_orient=None,
        ####SPI stuff ( ab initio)
        num_rings = 20,
        pixels_per_ring = 20,
        # first image orientation cardan angles [degrees] 
        SPI_x_rotation = 0,
        SPI_y_rotation = 0,
        SPI_z_rotation = 0,
        #crystallographic orientations (not consistent with SPI yet)
        # [ax_x,ax_y,ax_z] = vector parallel to rotation axis. Overridden if random orientations.        
        #num_orients_crys=100, # Miller indices orientations
        num_orients_crys=1, # Miller indices orientations
        #orientation_axis_crys = None,
        orientation_axis_crys = [1,1,1],
        #orientation_axis_crys = [0,0,1],#None,#[1,1,0]
        
        # for debugging/comparison with other works
        custom_cell_dims_for_miller_indices = None, #[17.174,14.93,13.384], # None # Implemented for comparison with others that use supercells.  
        override_max_q = False # False # Also special, implemented for comparison purposes but should be left as False by default.
        ######
    )
    same_deviations = False # whether same position deviations between damaged and undamaged crystal 

    # Optional: Choose previous folder for crystal results
    chosen_root_handle = None # None for new. use e.g. "tetra_v1", if want to add images under same params to same results.
    #=========================-------------------------===========================#

    ## DEBUG
    # WARNING we often assume that first crystal is damaged and second is undamaged when plotting. 
    crystal1_is_damaged = True # True  
    crystal2_is_damaged = False  # False



    #---------------------------Result handle names---------------------------#
    exp1_qualifier = "real"
    exp2_qualifier = "ideal"
    if chosen_root_handle is None:
        version_number = 1
        count = 0
        if tag != "":
            tag = "_" + tag        
        while True:
            if count > 299:
                raise Exception("could not find valid file in " + str(count) + " loops")
            results_parent_folder = RESULTS_LOCAL_PATH # needs to be synced with other functions
            root_handle = str(target) + tag
            exp_name1 = root_handle + "_" + exp1_qualifier + "_v" + str(version_number)
            exp_name2 = root_handle + "_" + exp2_qualifier  + "_v" + str(version_number)
            if path.exists(path.dirname(results_parent_folder + exp_name1 + "/")) or path.exists(path.dirname(results_parent_folder + exp_name2 + "/")):
                version_number+=1
                count+=1
                continue 
            break
    else:
        exp_name1 = chosen_root_handle + "_" + exp1_qualifier
        exp_name2 = chosen_root_handle + "_" + exp2_qualifier

    #exp_name2 = None

    #---------------------------------#
    water_index = None # None TODO automate
    pdb_path2 = None
    CNO_to_N = False
    S_to_N = False
    if target in["lys_salt","lys_no_salt","lys_salt_HF","lys_no_salt_HF","galliHigh", "lys_nass_probe_35"]:
        QUICK_TEST = False
        IDEAL = False
        COMPARE_REFINED = False
        include_H=False
        
        SPI = False
        #num_bragg_sets = 25
        #num_unique_supercells = 20
        #random_waters=702
        ####
        cycles_per_bragg_set = 1
        num_bragg_sets = 1
        num_unique_supercells = 50
        positional_stdv = 0
        num_supercells = 29*33*37
        #num_unique_supercells = 1
        #positional_stdv = 0.05
        #num_supercells = 29*33*37
        supercell_scale = 1
        random_waters=None
        zero_bfactors=False
        t_fineness=20
        crystal_qwargs["include_symmetries"]=True
        ###
        if QUICK_TEST or IDEAL:
            num_bragg_sets = 1
            num_unique_supercells = 1
            if QUICK_TEST:
                crystal_qwargs["include_symmetries"]=False
            if IDEAL:
                positional_stdv=0
        if COMPARE_REFINED:
            positional_stdv=0
            random_waters=None
            num_supercells=1
        #unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_1.5.hkl"
        unique_hkl ="/home/speno/AC4DC/scripts/scattering/targets/unique_reflections/unique_reflections_lysozyme_2.0.hkl"
        exp_qwargs["miller_indices_override"] = read_hkl(unique_hkl)
        exp_qwargs["spot_fraction_per_orient"] = 1/cycles_per_bragg_set
        exp_qwargs["num_orients_crys"] = cycles_per_bragg_set*num_bragg_sets
        exp_qwargs["t_fineness"]=t_fineness
        crystal_qwargs["supercell_simulations"] = num_unique_supercells
        crystal_qwargs["num_supercells"] = num_supercells
        crystal_qwargs["supercell_scale"] = supercell_scale
        crystal_qwargs["random_waters"] = random_waters
        crystal_qwargs["zero_bfactors"] = zero_bfactors
        crystal_qwargs["positional_stdv"]=positional_stdv
        laser_firing_qwargs["SPI"]=SPI

        target_handle = dict(
            lys_salt = "lys_salt_solvated_fast_H_5",
            lys_no_salt =  "lys_solvated_fast_H_6",
            lys_salt_HF =  "lys_salt_fast_high_fluence_2",
            lys_no_salt_HF =  "lys_solvated_fast_high_fluence_2",
            galliHigh = "lys_galli_HF_23",
            lys_nass_probe_35 = "nass_probe_35_1"
        )[target]        
        # TODO make this more sytematic.
        probe_delay = 0
        if target=="nass_probe_35":
            probe_delay = 35
        start_time += probe_delay
        end_time += probe_delay

        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et8.pdb" 
        if include_H:
            pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et8H.pdb" 
        if COMPARE_REFINED:
            pdb_path2 = dict(
                lys_salt = "/home/speno/AC4DC/scripts/scattering/targets/salt_group_1.pdb",
                lys_no_salt = "/home/speno/AC4DC/scripts/scattering/targets/no_salt_group_1.pdb", 
            )[target]
            crystal1_is_damaged = False
        else:
            assert(crystal1_is_damaged)
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/lys_8_cell.xpdb"; water_index = 69632
        CNO_to_N = False
        S_to_N = False
        folder = ""
        #allowed_atoms = ["H","C","N","O","S","Na","Cl","Gd_fast"]
        #allowed_atoms = ["H","C","N","O","S","Na","Cl"]
        allowed_atoms = ["C","N","O","S","Gd_fast"]

        if not laser_firing_qwargs["SPI"]:
            pass
            #exp_name2 = None # Don't do the undamaged target
    elif target == "neutze": #T4 virus lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/2lzm.pdb"
        target_handle = "lys-1_2"  
        folder = "lys" # If sim output folders are nested within subdir of __Molecular
        allowed_atoms = ["N_fast","S_fast"]
        CNO_to_N = True
    elif target == "hen": # egg white lys
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/4et8H.pdb"
        # Solvated targets
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/sol_4et8_full_struct_asym.xpdb"; water_index = 1089
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/sol_4et8_full_struct_unit_cell.pdb"; water_index = 8705
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/lys_asym_water.xpdb"; water_index = 
        #pdb_path = "/home/speno/AC4DC/scripts/scattering/solvate_1.0/lys_8_cell.xpdb"; water_index = 69632      
        # target_handle = "lys_nass_2"
        # folder = "lys"
        #'''
        # Light + Heavy atoms handles.
        target_handle = "lys_nass_Gd_full_1"#"lys_nass_Gd_full_1"#"lys_nass_no_S_3" #"lys_all_light-typical"#"lys_full-94_1"#"lys_nass_15"#"lys-5_3"#" #12keV, 0.1/0.01 count, 10 fs
        # Light atoms only
        #target_handle = "lys_nass_no_S_3"
        #folder = "lys" 
        folder = "" 
        #background_targets = "lys_water"
        #''' 
        '''
        target_handle = "lys_no_S_1"#"lys_no_S_2" #12keV, 0.1/0.01 count, 10 fs
        folder = ""        
        '''
        #//
        allowed_atoms = ["C","N","O"]
        #allowed_atoms = ["C","N","O"]
        #allowed_atoms = ["N","S_fast"]
        #allowed_atoms = ["N_fast"]
        #allowed_atoms = ["S_fast"]
        CNO_to_N = False
        S_to_N = True
        #//
    elif target == "tetra": 
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/5zck.pdb" 
        folder = ""#"tetra_CNO"
        target_handle = "lys_solvated_fast_high_fluence_2"#"lys_all_light-typical"#"6-5-2_tetra_CNO_3"
        #allowed_atoms = ["N_fast"]
        allowed_atoms = ["C","N","O"]
        CNO_to_N = False
        S_to_N = False
    elif target == "glycine":
        exp_qwargs["custom_cell_dims_for_miller_indices"] = [17.174,14.93,13.384]# # None,
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/glycine.pdb" 
        folder = ""
        target_handle = "lys_salt_fast_high_fluence_2" #"glycine_abdullah_high_H_6" #"lys_solvated_fast_high_fluence_2" # "glycine_abdullah_high_H_6" #"glycine_abdullah_4"
        allowed_atoms = ["C","N","O"]
        CNO_to_N = False
        S_to_N = False
    elif target == "copper_sulfate":
        allowed_atoms = ["Cu","S","O","H"]
        #allowed_atoms = ["Cu"]
        pdb_path = "/home/speno/AC4DC/scripts/scattering/targets/CuSO4.pdb" 
        target_handle = "copper_sulfate_above_e12_1fs_1" #"copper_sulfate_above_e12_14" #"copper_sulfate_below_e13_3#"copper_sulfate_above_e12_long_1"#"copper_sulfate_above_e12_14"
        folder = ""
        crystal_qwargs["cell_packing"]="triclinic"
        exp_qwargs["t_fineness"]=1
        slice_time = 1
        exp_qwargs["start_time"] = slice_time-0.1
        exp_qwargs["end_time"] = slice_time+0.1

    else:
        raise Exception("'target' invalid")
    #-------------------------------#
    if SEEDED:
        np.random.seed(0)
        
    import inspect
    src_file_path = inspect.getfile(lambda: None)
    sim_data_dir = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+folder)) + "/"


    # Set up experiments
    experiment1 = XFEL(exp_name1,energy,**exp_qwargs)
    experiment2 = XFEL(exp_name2,energy,**exp_qwargs)
    # Create Crystals

    crystal = Crystal(pdb_path,allowed_atoms,is_damaged=crystal1_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
    # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
    if pdb_path2 is None:
        if same_deviations:
            assert crystal_2_has_deviations
            # we copy the other crystal so that it has the same deviations in coords
            crystal2 = copy.deepcopy(crystal)#Crystal(pdb_path,allowed_atoms,cell_dim,is_damaged=False,CNO_to_N = CNO_to_N, **crystal_qwargs)
            crystal2.is_damaged = crystal2_is_damaged
        else:
            crystal2 = Crystal(pdb_path,allowed_atoms,is_damaged=crystal2_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
        if not crystal_2_has_deviations:
            crystal2.disable_pos_deviations()
        if show_crystal:
            crystal.plot_me(300000,water_index = water_index,template="plotly_dark")

    else:
        print("Setting second experiment to DIFFERENT crystal")
        crystal2 = Crystal(pdb_path2,allowed_atoms,is_damaged=crystal2_is_damaged,CNO_to_N = CNO_to_N,S_to_N=S_to_N, **crystal_qwargs)
#%
    if laser_firing_qwargs["SPI"]:
        SPI_result1 = experiment1.fire_laser(start_time,end_time,target_handle,sim_data_dir,crystal,results_parent_dir=results_parent_folder, **laser_firing_qwargs)
        SPI_result2 = experiment2.fire_laser(start_time,end_time,target_handle,sim_data_dir,crystal2,results_parent_dir=results_parent_folder,  **laser_firing_qwargs)
        stylin(exp_name1,exp_name2,experiment1.max_q,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2,custom_fig_width=fig_width,custom_fig_height=fig_height)
    else:
        experiment1.fire_laser(start_time,end_time,target_handle,sim_data_dir,crystal, results_parent_dir=results_parent_folder, **laser_firing_qwargs)
        exp1_orientations =experiment1.used_orientations
        create_reflection_file(exp_name1,results_parent_dir=results_parent_folder)
        rfl_to_sca(exp_name1)
        if exp_name2 != None:
            laser_firing_qwargs["random_orientation"] = False
            experiment2.set_orientation_set(exp1_orientations)  # pass in orientations to next sim, random_orientation must be false!
            experiment2.fire_laser(start_time,end_time,target_handle,sim_data_dir,crystal2, results_parent_dir=results_parent_folder, **laser_firing_qwargs)
            create_reflection_file(exp_name2,results_parent_dir=results_parent_folder)
            rfl_to_sca(exp_name2)

        stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,custom_fig_width=fig_width,custom_fig_height=fig_height) # Note we are passing the max q, not max q_scr.
#^^^^^^^
#%% Pixels/SPI
if __name__ == "__main__":
    #fig_width = 3.49751 # 20
    #fig_height = fig_width*3/4 # 20
    stylin(exp_name1,exp_name2,experiment1.max_q,SPI=laser_firing_qwargs["SPI"],SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2, custom_fig_height=fig_height,custom_fig_width=fig_width,
           min_R_dmg_pixel=0,spi_full_rings_only=False,log_range=None,
           log_diff_vmin = -0.4, 
           log_diff_vmax = 0.6, dpi=800,)
#%% Miller/macrocrystal
if interactive and __name__ == "__main__":
    stylin(exp_name1,exp_name2,experiment1.q_to_X(experiment1.max_q)/1e7,show_labels=False) # Note we are passing the max q, not max q_scr.
    #stylin("glycine_v36__real","glycine_v36__ideal",2.3)
    #stylin(exp_name1,exp_name2,3)
#%%--------STRUCTURE CONSTRUCTOR------

# Save full structures in pdb format for SOLVATE
# Using this structure is not amazing practice, it takes a lot of time and potentially memory!
# SOLVATE allows for generating just the water with the solute removed. So an alternative method might 
# be generating the water for an individual unit cell with different distributions (seems possible by using slightly different thickness), 
# stitching them together, and removing atoms that are outside the Wigner–Seitz cell.
# We then calculate the form factor for the crystal, followed by the form factor for each of N water cells by defining water_background = Crystal(water_background_N,allowed_atoms).
# Finally, the rest of the water drop could be calculated by generating a large distribution of water, then scaling its contribution to the form factor.
if interactive and __name__ == "__main__":
    ##### Crystal params
    #pdb_file ="4et8H.pdb"
    #pdb_file ="2qspH.pdb"
    pdb_file ="9EPD_Hfix_singleconf.pdb"
    targets_dir = path.abspath(path.join(__file__ ,"../")) + "/targets/"
    pdb_path = targets_dir + pdb_file
    crystal_qwargs = dict(
        supercell_scale = 1,  # for SC: supercell_scale^3 unit cells
        positional_stdv = 0, 
        include_symmetries = True,  # should unit cell contain symmetries or just one asymmetric unit?
        cell_packing = "triclinic",#"SC",
    )
    custom_residue_name=None
    #allowed_atoms = ["C","N","O","S"]
     #I3C
    # custom_residue_name="I3C"
    #allowed_atoms = ["C","N","O","I","H"]
    #allowed_atoms =["O","S","H","Cu"]
    #allowed_atoms =["C","N","O","H","S","Na","Cl","Gd"]
    allowed_atoms=get_sim_elements("hemoglobin_2QSP-9")
    crystal = Crystal(pdb_path,allowed_atoms,is_damaged=False, **crystal_qwargs)
    #crystal.save_structure(custom_residue_name=custom_residue_name,chain_name="A")
    crystal.save_structure_by_reference()

    
   

# %%
def plot_recovered_atoms():
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
if interactive and __name__ == "__main__":
    plot_recovered_atoms()    
# %%



