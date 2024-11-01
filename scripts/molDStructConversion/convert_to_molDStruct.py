#%%
import sys
import os
import pandas as pd
import os.path as path
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
from scatter import XFEL,Crystal,stylin
from core_functions import get_sim_params,get_pdb_path
import imaging_params as imaging_params
import numpy as np

# Converts AC4DC data to IONIZATION_DATA used for input to MolDStruct CR-MD.
# Also creates pdb file.

# Note using this requires a modified version of the reading of data in md.c

target_path = ""


crystal_params = dict(
    supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
    num_supercells = 1,
    supercell_simulations = 1,        
    include_symmetries = True, ##True  # should unit cell contain symmetries?
    positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
    cell_packing = "SC",
    rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
)


def get_random_charge_states(element):
    SEEDED = False
    element.times_used = element.crystal.ff_calculator.get_times_used()
    if element.num_atoms != len(element.crystal.sym_rotations)*len(element.coords):
        raise Exception("num atoms was not same on set_stochastic_states call as when set by set_coord_deviation")
    if element.crystal.is_damaged:
        charges = np.empty(shape = (element.num_atoms,len(element.times_used)))  # (num atoms, times)   
        for idx in range(element.num_atoms):
            seed = None
            if SEEDED:
                seed = idx
            charges[idx] = element.crystal.ff_calculator.random_charge_snapshots(element.name,seed) 
    return charges



def create_charge_file(charges,element,subdir_name,output_dir,overwrite=False):
    '''
    Generates a
    '''
    print(f"Creating {element} charge file for {subdir_name}")
    df = pd.DataFrame(charges)
    # Save as file
    save_dir = output_dir+subdir_name+"/"
    os.makedirs(save_dir, exist_ok=True) 
    save_path = save_dir+element+'_charges.csv'
    if os.path.isfile(save_dir): 
        if not overwrite:
            print("Cannot write, file already present at",save_dir)
            return
        os.remove(save_dir)

    df.to_csv(save_path, index=False)



#PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"I3C") 

SCRIPTS_PATH = path.abspath(path.join(__file__ ,"../../"))

SCATTER_DIR = SCRIPTS_PATH +"/scattering/"
PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"tetra") 
OUTPUT_PATH =  path.abspath(path.join(__file__ ,"../"))+ "/output/"

MOLECULAR_PATH = path.abspath(path.join(SCRIPTS_PATH, "../output/__Molecular/")) + "/" # directory of damage sim output folders


allowed_atoms = ["C","N","O","I"]



sim_handle = "nass_probe_37_2"

crystal = Crystal(PDB_STRUCTURE,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True,**crystal_params)

param_dict,_,_ = get_sim_params(sim_handle)
start_time = param_dict["start_t"]
end_time = param_dict["end_t"]
energy = param_dict["energy"]



# Assign plotter object to calculate charges (because code debt)
xfel = XFEL("dummy",energy)
ff_calculator = xfel.get_ff_calculator(start_time,end_time,sim_handle,MOLECULAR_PATH)     
crystal.set_ff_calculator(ff_calculator)    


times_used = None

species_charges = {}
for element in crystal.species_dict.keys():
    species_charges[element] = get_random_charge_states(crystal.species_dict[element])


# Generate csv file
for element, charges in species_charges.items():
    PDB_element = element.split("_")[0]
    create_charge_file(charges,PDB_element,"test",OUTPUT_PATH)




# class Target(Enum):
#     UNIT = imaging_params.goldilocks_dict_unit
#     NINE = imaging_params.goldilocks_dict_3x3x3


#im_params = Target.NINE
#crystal = Crystal(PDB_STRUCTURE,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **im_params["crystal"])

# %%
