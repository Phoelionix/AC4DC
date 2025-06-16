# Because I didn't implement this when loading backups. (TODO)
#%%
import stringprep
import matplotlib.rcsetup as rcsetup
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import matplotlib.colors as colors
import sys
import glob
import csv
import subprocess
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress
from core_functions import get_mol_file, parse_elecs_from_latex, get_sim_params, ATOMS, ATOMNO
from scipy.interpolate import splrep, splev
from scipy.signal import savgol_filter
import pandas as pd
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)

def StitchTogetherDeltaOccupancies(handles,out_handle,molecular_path = None):
    if molecular_path is None:
        molecular_path = path.abspath(path.join(__file__ ,"../../output/__Molecular/")) + "/"

    # MUST be in order

    outDir = molecular_path+out_handle
    os.makedirs(outDir,exist_ok=True)
    delta_files_dict = {}
    for handle in handles:
        handleDir = molecular_path + handle
        delta_files = glob.glob(f"{handleDir}/delta_occupancy_from_*.csv")
        for delta_file in delta_files:
            K = os.path.basename(delta_file)
            if K not in delta_files_dict:
                delta_files_dict[K] = [delta_file]
            else:
                delta_files_dict[K].append(delta_file)


    for delta_file_basename, delta_file_list in delta_files_dict.items():
        out_path = f"{outDir}/{delta_file_basename}"
        with open(out_path,'w') as f:
            f.write("# Ionic changes (stitched)\n# Time (fs) | transition to charge")
        
        raw_data_list = []
        for delta_file in delta_file_list:
            with open(delta_file,'r') as f:
                lines = f.readlines()
                raw_data_list.append(np.genfromtxt(lines, comments='#', dtype=np.float64))
        max_time = None
        combined_data = None
        for i, raw_data in enumerate(raw_data_list):
            first_time_to_be_found = True
            for t_idx, line in enumerate(raw_data):
                if np.sum(line[1:])>0:
                    if first_time_to_be_found:
                        min_time_idx = t_idx
                        assert max_time is None or max_time < line[0]
                    first_time_to_be_found = False
                    max_time=line[0]
                    max_time_idx = t_idx
            if combined_data is None:
                combined_data = raw_data[:max_time_idx+1]
            else: 
                np.append(combined_data,raw_data[min_time_idx:max_time_idx+1])
        
        # the last file is the times we want
        target_times = raw_data_list[-1][:, 0] - 1e-10
        idxes = np.searchsorted(combined_data[:,0],target_times)
        idxes[idxes==len(combined_data)]-=1
        combined_data = combined_data[idxes,:]
            

        with open(out_path, "ab") as f:
            str_data = np.where(
                combined_data==0,
                '0',
                np.char.mod('%.5e', combined_data)
            )
            np.savetxt(f, str_data, fmt='%s')


if __name__ == "__main__":
    handles = ["I3C_25fs_backup","I3C_25fs_backup2","I3C_25fs_1"]
    out_handle = "I3C_stitched"
    StitchTogetherDeltaOccupancies(handles,out_handle)


            




    # No duplicates!
# %%
