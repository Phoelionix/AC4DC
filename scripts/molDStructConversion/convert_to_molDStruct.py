#%%
import sys
import os
import pandas as pd
import os.path as path
sys.path.append('/home/speno/AC4DC/scripts/pdb_parser')
sys.path.append('/home/speno/AC4DC/scripts/scattering')
sys.path.append('/home/speno/AC4DC/scripts/')
from scattering.scatter import XFEL,Crystal,stylin,Atomic_Species
from core_functions import get_sim_params,get_sim_elements,get_pdb_path,ATOMNO
import scattering.imaging_params as imaging_params
import numpy as np
from scipy import constants as C
import struct
import textwrap
from plotter_core import Plotter # For typing
import matplotlib.pyplot as plt


# Converts AC4DC data to IONIZATION_DATA used for input to MolDStruct CR-MD.

# sim_handle = "tmp_I3C_2"
# num_steps = 500
#sim_handle = "I3C_55fs_4"
#num_steps = 4900   # NOT in attoseconds when doing i3c 55 fs
#target = "I3C.gro"
#sim_handle = "lys_salt_solvated_fast_H_4"
AVERAGE_CHARGES = None 
ALLOW_SELECT_SAME_TIMES = True
SAVE_CSV_COPY = True


def get_charge_states(element:Atomic_Species,element_charge_snapshot_selector=None):
    if element_charge_snapshot_selector is None:
        #element_charge_snapshot_selector = element.crystal.ff_calculator.random_charge_snapshots
        element_charge_snapshot_selector = element.crystal.ff_calculator.continuity_charge_snapshots
    SEEDED = False
    element.times_used = element.crystal.ff_calculator.get_times_SCATTER()
    if element.get_num_atoms() != len(element.crystal.sym_rotations)*len(element.coords):
        raise Exception("num atoms was not same on set_stochastic_states call as when set by set_coord_deviation")
    if element.crystal.is_damaged:
        charges = np.empty(shape = (element.get_num_atoms(),len(element.times_used)),dtype=int)  # (num atoms, times)   
        for idx in range(element.get_num_atoms()):
            seed = None
            if SEEDED:
                seed = idx
            charges[idx] = element_charge_snapshot_selector(element.name,seed) 
            #print(charges[idx])
            
            # TEMPORARY HACK COS UNSIGNED SHORT DUMBNESS
            # (Can't pass negative values...)
            #if element.name == "I_fast":
            charges[idx] = np.maximum(0,charges[idx])
            if not np.all(charges.astype(np.ushort)[idx] <= ATOMNO[element.name]):
                print(charges.astype(np.ushort)[idx])
                print(ATOMNO[element.name])
                raise Exception("charge greater than atomic number!")

    return charges

def binary(num):
    return ''.join('{:0>8b}'.format(c) for c in struct.pack('!f', num))

def create_charge_file(charges,element,save_dir,overwrite=False,csv=False):
    '''
    Generates a

    '''
    save_dir_csv = save_dir
    save_dir_bin = save_dir + "IONIZATION_DATA/"
    print(f"Creating {'full' if element is None else element} charge file in {save_dir}")
    os.makedirs(save_dir_bin, exist_ok=True) 
    os.makedirs(save_dir_csv, exist_ok=True) 
    if element is None:
        save_path = f"{save_dir_bin}charges"
        save_path_csv = f"{save_dir_csv}charges"
    else:
        save_path = f"{save_dir_bin}{element}_charges"
        save_path_csv = f"{save_dir}{element}_charges"
    
    if os.path.isfile(save_dir): 
        if not overwrite:
            print("Cannot write, file already present at",save_dir)
            return
        os.remove(save_dir)

    # with open(save_path+".bin", "wb") as file:
    #     newFileByteArray = bytearray(charges)
    #     file.write(newFileByteArray)
    charges.astype(np.ushort).swapaxes(0,1).tofile(save_path+".bin") #uint16, shape = (num timesteps, num atoms)
    
    if csv:
        df = pd.DataFrame(charges)
        df.to_csv(save_path_csv+".csv", index=False)

    with open(save_path+".bin", "rb") as file:
        n=10
        print(f"First {n} binary data elements:")
        print(struct.unpack('H'*n, file.read(2*n)))


def create_data_file(data,tag,save_dir,overwrite=False,csv=False):
    '''
    Generates a
    '''
    save_dir_csv = save_dir
    save_dir_bin = save_dir + "IONIZATION_DATA/"
    
    data[0] = data[1]  # Patch.

    print(f"Creating {tag} file in {save_dir.split('/')[-1]}")
    os.makedirs(save_dir_bin, exist_ok=True) 
    os.makedirs(save_dir_csv, exist_ok=True) 
    save_path = save_dir+"IONIZATION_DATA/"+tag
    save_path_csv = save_dir+tag
    if os.path.isfile(save_dir): 
        if not overwrite:
            print("Cannot write, file already present at",save_dir)
            return
        os.remove(save_dir)

    # with open(save_path+".bin", "wb") as file:
    #     newFileByteArray = bytearray(data)
    #     file.write(newFileByteArray)
    data.astype('float32').tofile(save_path+".bin")
    
    if csv:
        df = pd.DataFrame(data)
        df.to_csv(save_path_csv+".csv", index=False)

    with open(save_path+".bin", "rb") as file:
        n=10
        print(f"First {n} binary data elements:")
        print(struct.unpack('f'*n, file.read(4*n)))

#PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"I3C") 

SCRIPTS_DIR = path.abspath(path.join(__file__ ,"../../")) + "/"

SCATTER_DIR = SCRIPTS_DIR +"scattering/"
#PDB_STRUCTURE = get_pdb_path(SCATTER_DIR,"tetra") 
OUTPUT_PATH =  path.abspath(path.join(__file__ ,"../"))+ "/output/"

TARGET_DIR = SCATTER_DIR+ "targets/"


# already defined in plotter_core
#MOLECULAR_PATH = path.abspath(path.join(SCRIPTS_DIR, "../output/__Molecular/")) + "/" # directory of damage sim output folders


# def LennardJones():
#     out_folder = get_save_folder()
#     f = open(out_folder+"/lennard_jones_parameters.txt", "w") 

#     for i, j in zip(types, atom_number):
#         f.writelines("{0} 0 0\n".format(j))
#     f.close()

RANDOM_CHARGES = True # TEMPORARY

def get_save_folder(sim_handle):
    if AVERAGE_CHARGES:
        tag = "avg"
    elif RANDOM_CHARGES:
        tag = "stoch"
    else: 
        tag = "stoch_cont"
    folder_name = f"{sim_handle}_{tag}"
    return OUTPUT_PATH + folder_name + "/"



def charges(crystal:Crystal,ff_calculator:Plotter,sim_handle,csv=False,individual_elements = False,average_charges=AVERAGE_CHARGES):
    if average_charges is None:
        average_charges = True 
    print("Beginning writing of charges...")

    out_folder = get_save_folder(sim_handle)
    print(OUTPUT_PATH)
    
    num_steps = len(ff_calculator.get_times_SCATTER())
    species_charges = {}
    num_atoms = 0
    num_atoms_for_print = 0
    print(f"Processing charges at {num_steps} time steps for:")
    for element in crystal.species_dict.keys():
        num_atoms_for_print+=len(crystal.species_dict[element].serial_numbers)
        print(f"{len(crystal.species_dict[element].serial_numbers)} {element} atoms")
    print(f"Total: {num_atoms_for_print}")

    print("Generating:")
    for element in crystal.species_dict.keys():
    #for element in ["Na","Cl","S","H","C","N","O"]:
        print(f"{element}...")
        element_obj:Atomic_Species = crystal.species_dict[element]
        if average_charges:
            species_charges[element] = element_obj.crystal.ff_calculator.get_average_charge_ff_calculator(element_obj.name)
        else:
            charge_selector = None
            if RANDOM_CHARGES:
                charge_selector = element_obj.crystal.ff_calculator.random_charge_snapshots
            species_charges[element] = get_charge_states(element_obj,element_charge_snapshot_selector=charge_selector)
            # TEMPORARY HACK COS UNSIGNED SHORT DUMBNESS
            # (Can't pass negative values...)
            #if element_obj.name == "I_fast":
            species_charges[element] = np.maximum(0,species_charges[element])
        num_atoms += element_obj.get_num_atoms()
    assert num_atoms > 0
    
 
    # Order charges in order that matches the structure file. 
    combined_charges = np.empty(shape = (num_atoms,num_steps))  # (num atoms, times)   
    species_list = np.empty(shape = (num_atoms,),dtype=object)
    tuples=[]
        
    # Need to track the serial numbers so that we can create charge file in order, while also allowing for possibility that there are jumps in serial num. 
    j=0
    for element, charges in species_charges.items():
        for i, s_num in enumerate(crystal.species_dict[element].serial_numbers):
            if average_charges:
                tuples.append((charges,element,s_num))
            else:
                tuples.append((charges[i],element,s_num))
            j+=1
        tuples.sort(key=lambda x: x[2]) # sort by serial num
        for i, (charges,_,_) in enumerate(tuples):
            combined_charges[i]=charges


        if individual_elements:
            PDB_element = element.split("_")[0]
            create_charge_file(charges,PDB_element,out_folder,csv=csv)    # Each row is an atom. each column is a time step.

        # for q in charges:
        #     assert q < ATOMNO[element], f"{q},{ATOMNO[a]}"
        
    create_charge_file(combined_charges,None,out_folder,csv=csv)    # Each row is an atom. each column is a time step.
    for charges, a, _ in tuples:
        for q in charges:
            assert q <= ATOMNO[a], f"{q},{ATOMNO[a]}"


def DebyeLength(ff_calculator:Plotter,sim_handle,csv=False):
    pl = ff_calculator
    

    tempList = []
    denseList = []
    for t in pl.get_times_SCATTER():
        tempList.append( pl.get_temp(t, 1000) ) # eV
        # denseList.append( pl.get_free_electron_density(t) ) # per angstrom cube

    T = np.array(tempList)
    n = T [: ,1]
    T = T[:,0]
    print("n:",n)
    print("T:",T)

    # Can ignore if at step 0
    for i, elem in enumerate(T):
        if elem <= 0:
            print(f"Warning, T at step {i} = {elem}")
    for i, elem in enumerate(n):
        if elem <= 0:
            print(f"Warning, n at step {i} = {elem}")

    eps = 1e-12
    n = np.maximum(eps,n)
    lambdaD=np.sqrt(C.epsilon_0 * C.nano * T *C.eV / n /C.e/C.e) # should have units nm
    #lambdaD=np.sqrt(C.epsilon_0 * C.angstrom * T *C.eV / n /C.e/C.e) # should have units Angstrom

    out_folder = get_save_folder(sim_handle)
    create_data_file(T*11606,"electron_temperature",out_folder,csv=csv) # K
    create_data_file(n/C.nano**3,"electron_density",out_folder,csv=csv) # nm^-3
    create_data_file(lambdaD,"debye_data",out_folder,csv=csv) # nm




def get_plotter(handle,parent_dir_path,start_time,end_time,t_fineness)->Plotter:
    ff_calculator = Plotter(handle,parent_dir_path,out_prefix_text = "Setting up plotter...",
                            skip_mol_file=True,
                            initialise=False)
    plt.close()
    ff_calculator.initialise_charges_only()
    ff_calculator.initialise_form_factor_params(start_time,end_time,None,None,t_fineness=t_fineness) # q_fineness isn't used for our purposes.   
    return ff_calculator

def convert_to_molDStruct(sim_handle:str,target_path:str,num_steps:int,allowed_atoms,start_time,end_time,debye=True,sim_parent_dir_path=None,param_dict=None):


    crystal_params = dict(
        supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
        num_supercells = 1,
        supercell_simulations = 1,        
        include_symmetries = None, ##True  # should unit cell contain symmetries?
        positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
        cell_packing = "SC",
        rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
    )

    crystal = Crystal(target_path,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=False,**crystal_params)
    # Assign plotter object to calculate charges (this is code debt)
    ff_calculator = get_plotter(sim_handle,sim_parent_dir_path,start_time,end_time,t_fineness=num_steps)   
    ff_calculator.allow_select_same_times = ALLOW_SELECT_SAME_TIMES  
    crystal.set_ff_calculator(ff_calculator) 


    os.makedirs(get_save_folder(sim_handle), exist_ok=True) 
    
    charges(crystal,ff_calculator,sim_handle,csv=SAVE_CSV_COPY)
    if debye:
        DebyeLength(ff_calculator,sim_handle,csv=SAVE_CSV_COPY)


    log_file = get_save_folder(sim_handle) + "log.txt"
    target_handle = '.'.join(os.path.basename(target_path).split('.')[:-1])
    with open(log_file,'w') as f:
        f.write(textwrap.dedent(f"""\
                target: {target_handle}
                plasma simulation: {sim_handle}
                p. sim. parameters: {param_dict if param_dict is not None else '-'}
                """))
    
    print("Done! Remember to sit straight!")

def convert_to_molDStruct_standard(sim_handle:str,target_path:str,num_steps:int,debye=True,sim_parent_dir_path=None):
        allowed_atoms = get_sim_elements(sim_handle,molecular_path=sim_parent_dir_path)


        param_dict,_,_ = get_sim_params(sim_handle,molecular_path=sim_parent_dir_path)
        start_time = param_dict["start_t"]
        end_time = param_dict["end_t"]
        #energy = param_dict["energy"]




        convert_to_molDStruct(sim_handle,target_path,num_steps,allowed_atoms,start_time,end_time,debye,sim_parent_dir_path,param_dict)


# def convert_to_molDStruct_standard(sim_handle:str,target_path:str,num_steps:int,debye=True,sim_parent_dir_path=None):
#         allowed_atoms = get_sim_elements(sim_handle,molecular_path=sim_parent_dir_path)


#         crystal_params = dict(
#             supercell_scale = 1,  ##3 # for SC: cell_scale^3 unit cells 
#             num_supercells = 1,
#             supercell_simulations = 1,        
#             include_symmetries = None, ##True  # should unit cell contain symmetries?
#             positional_stdv = 0, # Introduces disorder to positions. Note this is a deviation from the IDEAL structure, so is not a measure of similarity with undamaged and damaged structure but the ideal structure to recover and the dmaaged structure. Can roughly model atomic vibrations/crystal imperfections. Should probably set to 0 if quickly gauging serial crystallography R factor, as should somewhat average out.
#             cell_packing = "SC",
#             rocking_angle = 1.2,  # (approximating mosaicity, infinite crystal sim only)
#         )

#         crystal = Crystal(target_path,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=False,**crystal_params)


#         target_handle = '.'.join(os.path.basename(target_path).split('.')[:-1])

#         param_dict,_,_ = get_sim_params(sim_handle,molecular_path=sim_parent_dir_path)
#         start_time = param_dict["start_t"]
#         end_time = param_dict["end_t"]
#         energy = param_dict["energy"]



#         # Assign plotter object to calculate charges (because code debt)
#         xfel = XFEL("dummy",energy,t_fineness=num_steps)
#         ff_calculator = xfel.get_ff_calculator(start_time,end_time,sim_handle,sim_parent_dir_path)   
#         ff_calculator.allow_select_same_times = ALLOW_SELECT_SAME_TIMES  
#         #crystal.plot_me()
#         crystal.set_ff_calculator(ff_calculator)    

#         # TODO json format.
#         log_file = get_save_folder(sim_handle) + "log.txt"
#         os.makedirs(get_save_folder(sim_handle), exist_ok=True) 

#         charges(crystal,ff_calculator,csv=SAVE_CSV_COPY)
#         if debye:
#             DebyeLength(ff_calculator,csv=SAVE_CSV_COPY)

#         with open(log_file,'w') as f:

if __name__ == "__main__":

    #sim_handles = ["lys_salt_solvated_fast_H_4","lys_solvated_fast_H_4"]
   # sim_handles = ["lys_solvated_H_2",]
    #sim_handles = ["lys_solvated_H_9","lys_salt_solvated_H_1"]
    #sim_handles = ["I3C_25fs_combined"]
    #sim_handles = ["lys_salt_fast_high_fluence_2"]

    #num_steps = 3600 #4900 #3600 #
    if len(sys.argv)!=3:
        print("Usage: python3.9 sim_output_handle path/to/gro/file")
        quit()
    sim_handle,target_path=sys.argv[1:3]

    sim_params=get_sim_params(sim_handle)[0]
    sim_duration_fs=sim_params["end_t"]-sim_params["start_t"]

    num_steps = int(np.ceil(50*sim_duration_fs))
    print(f"Creating charge file with {num_steps} steps")

    #target = "CNO_debug.gro"
    #target = "4et8.gro"
    #target = "lys_example.gro"
    #target = "4et8H_full_struct_Hfix.gro"
    #target_path=TARGET_DIR + target
    AVERAGE_CHARGES = False

    convert_to_molDStruct_standard(sim_handle,target_path,num_steps)






    # class Target(Enum):
    #     UNIT = imaging_params.goldilocks_dict_unit
    #     NINE = imaging_params.goldilocks_dict_3x3x3


    #im_params = Target.NINE
    #crystal = Crystal(PDB_STRUCTURE,allowed_atoms,is_damaged=True,convert_excluded_elements_to_N=True, **im_params["crystal"])

# %%






