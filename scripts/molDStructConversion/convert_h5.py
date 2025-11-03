#%%
import h5py
import numpy as np
import os,sys
import os.path as path
import pandas as pd
from convert_to_molDStruct import create_data_file,convert_to_molDStruct,get_save_folder


#####
ATOMS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'
       +' Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe').split()
ATOMNO = {}
i = 1
for symbol in ATOMS:
    ATOMNO[symbol] = i
    i += 1
#####

SCRIPTS_DIR = path.abspath(path.join(__file__ ,"../../")) + "/"

SCATTER_DIR = SCRIPTS_DIR +"scattering/"
OUTPUT_PATH =  path.abspath(path.join(__file__ ,"../"))+ "/output/"

TARGET_DIR = SCATTER_DIR+ "targets/"
TARGET_PATH= TARGET_DIR+"hemoglobin_solv_Hfix.gro"


def data_to_ac4dc_format(times,values):
    # AC4DC format: time is first element of each row.
    assert times.shape==(values.shape[0],)
    #shaped_times=np.reshape(times,(times.shape[0],)+(1,)*(len(values.shape)-1))

    return np.column_stack((times,values))

def gen_files(path,runid,mimic_sim_out_dir):
    with h5py.File(path, "r") as f:
        runname = f"run{runid:02d}"
        #runname = "run%02d" % runid
        fluence = f["fluence"][()][runid]
        stime = 1e15 * f["stime"][()]
        free_electron_temperature = f[runname]["tev"][()]
        free_electron_density = f[runname]["ne"][()]
        debye_length = f[runname]["dl"][()]

        n_steps = stime.shape[0]  

        #df = pd.DataFrame(stime)
        
        charge_states_dict={}
        elements_dict = {1:"H", 
                    2:"C",
                    3:"N",
                    4:"O",
                    5:"Na",
                    6:"S",
                    7:"Cl",
                    8:"Fe"}
        
        for key,ele in elements_dict.items():
            print(f"saving {key,ele}" )
            charge_states = f[runname][f"yiso_{key}"]
            assert charge_states.shape[-1]==ATOMNO[ele]+1

            # $1s^{2}2p^{8}3p^{16}4p^{2}5p^{3}$ $1s^{2}2p^{8}3p^{16}4p^{2}5p^{2}$
            orbitals_header = "#           | " + " ".join(["1s^{"+f"{ATOMNO[ele]-n}"+"}$" for n in range(ATOMNO[ele]+1)] )
            print(orbitals_header)

            charge_states=data_to_ac4dc_format(stime,charge_states)
            mimic_handle=f"converted_charges_{os.path.basename(path).split('.')[0]}_{runid}"
            converted_data_folder=f"{mimic_sim_out_dir}/{mimic_handle}/"
            os.makedirs(converted_data_folder,exist_ok=True)
            np.savetxt(f"{converted_data_folder}/dist_{ele}.csv", charge_states, delimiter=" ",header=orbitals_header)
            charge_states_dict[ele]=charge_states

        start_time,end_time= stime[0],stime[-1]
        convert_to_molDStruct(sim_handle=mimic_handle,target_path=TARGET_PATH,num_steps=n_steps,
            allowed_atoms=elements_dict.values(),start_time=start_time,end_time=end_time,
            debye=False,sim_parent_dir_path=mimic_sim_out_dir
            ) 

        # do Debye separately since already have the data
        out_dir=get_save_folder(mimic_handle)
        create_data_file(free_electron_temperature*11606,"electron_temperature",out_dir,csv=False) # eV --> K
        create_data_file(1e-9*free_electron_density,"electron_density",out_dir,csv=False) # cm^-3 --> nm^-3
        create_data_file(1e9*debye_length,"debye_data",out_dir,csv=False) # m --> nm
            

            # df = pd.DataFrame(charge_states)
            # df.to_csv("path/to/file.csv",header=False,index=False)
            
#for handle in ["3fs","10fs"]:
# for handle in ["10fs"]:
#     gen_files(f"hdf5_files/{handle}.h5",f"{handle}_converted/")

if __name__ == "__main__":
    assert len(sys.argv)==4
    hdf5_folder_path,handle,runid= sys.argv[1:]
    gen_files(os.path.join(hdf5_folder_path,f"{handle}.h5"),int(runid),
              mimic_sim_out_dir=os.path.join(path.abspath(path.join(__file__ ,"../")),"output",f"{handle}_{runid}_converted","")) 

    #np.savetxt('output.csv', data, delimiter=',', fmt='%f')
# %%
