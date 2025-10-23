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
from xpdb import sloppyparser as xPDBParser
from xpdb import SloppyPDBIO as xPDBIO
from xpdb import SloppyStructureBuilder as xStructureBuilder  # Hack, enables atom counts over 10,000
from Bio.PDB.Atom import Atom as PDB_Atom
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
import colorcet as cc; import cmasher as cmr
from mpl_toolkits.mplot3d import Axes3D
import plotly.graph_objects as go
import plotly.offline as pltly_offline
from IPython.display import display, HTML
from IPython import get_ipython
from core_functions import get_sim_elements
interactive = False
if interactive and __name__ == "__main__":
    get_ipython().run_line_magic('colors', 'nocolor')
    get_ipython().run_line_magic('matplotlib', 'widget')
    # pltly_offline.init_notebook_mode()
    # display(HTML(
    #     '<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.1/MathJax.js?config=TeX-MML-AM_SVG"></script>'
    # ))   

import pandas as pd
import csv

import glob
import random 

from scatter import rfl_to_sca, Results, ScalingByCopyingSigmaRatio

import inspect



def all_reflections_to_scalepack(result_handle,results_parent_dir,out_dir=None, reflections_dir=None,tag_override=None,scaling_method=None,scaling_thing=False,create_mtz=False):
    sample_results_and_create_scalepack(result_handle,"ALL",results_parent_dir,out_dir=out_dir,reflections_dir=reflections_dir, tag_override=tag_override,scaling_thing=scaling_thing,create_mtz=create_mtz)
    
def sample_results_and_create_scalepack(result_handle,num_to_sample,results_parent_dir,out_dir=None, reflections_dir=None,tag_override=None,scaling_method=None,scaling_thing=False,create_mtz=False):
    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    if out_dir is None:
        out_dir=scattering_dir+"random_sample/scalepack/"
    if reflections_dir is None:
        reflections_dir = "random_sample/reflections/"
        
    reflections_handle= sample_results_and_create_reflection_file(result_handle,num_to_sample,results_parent_dir,out_directory=reflections_dir,tag_override=tag_override)
    reflections_handle = reflections_handle.split("/")[-1]
    if scaling_thing:
        # TESTING TEMPORARY TODO
        scaling_method=None
        cif_file = f"/home/speno/PhenixWorkspace/data/4et8-sf.cif"
        scaling_method = ScalingByCopyingSigmaRatio(cif_file,scattering_dir+"random_sample/reflections/"+ reflections_handle + ".rfl")
        #
    rfl_to_sca(reflections_handle,reflections_dir,out_dir,scaling_method=scaling_method,create_mtz=create_mtz)
    rfl_to_sca(reflections_handle+"_unmerged",reflections_dir,out_dir,scaling_method=scaling_method,create_mtz=create_mtz)
    
    #reflections_dir=scattering_dir+"random_sample/reflections/"
    #out_directory=scattering_dir+"random_sample/scalepack/"

def sample_results_and_create_reflection_file(result_handle,num_to_sample,results_parent_dir,out_directory="random_sample/reflections/",tag_override=None):
    '''
    Generates a **merged** .rfl file by sampling all result files in directory given by result_handle, including subdirectories

    '''
    print("Creating reflection file for",result_handle)
    results_dir = results_parent_dir+ result_handle+"/"
    assert path.isdir(results_dir), f"Directory not found: {results_dir}" 
    os.makedirs(out_directory, exist_ok=True) 
    init = False

    filenames = glob.glob(results_dir+"/**/*.pickle",recursive=True)
    if num_to_sample == "ALL":
        num_to_sample = len(filenames)
    assert num_to_sample<=len(filenames), f"Only found {len(filenames)} files, but trying to sample {num_to_sample}"

    # Sample num_to_sample result files
    random.shuffle(filenames)
    for filename in filenames[0:num_to_sample]:
        result = Results.get_result(filename,"")[0]
        if result == "__PASS__":
            continue
        if result == None:
            continue    

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
    if tag_override is None:
        tag = 1
        out_path = f"{out_directory}{result_handle}_{tag}"
        # Iterate until get unique tag
        while glob.glob(f"{out_path}*"): 
            tag += 1
            out_path = f"{out_directory}{result_handle}_{tag}"
            if tag > 9999:
                raise Exception(f"Over {tag} folders with same handle in output directory")
    else: 
        out_path = f"{out_directory}{result_handle}"
        if tag_override !="":
            out_path +=f"_{tag_override}"
    #columns.reverse() #???
    df = df.sort_values(by=["l","k","h"],axis=0)
    for i in ("hkl"):
        df[i] = df[i].astype('int')
    df["I"] = df["I"].astype('float')
    df = df.round(6)
    #df.drop_duplicates(subset = ["h","k","l"],inplace=True) # TODO should take average.
    df = df[df['I']>=0.01] # TODO temporary fix for appearance of low values that needs to be squashed.
    df.to_csv(out_path+"_unmerged.rfl",header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")

    df_merged = df.groupby(["h","k","l"]).mean().reset_index()
    df_merged = df_merged.sort_values(by=["l","k","h"],axis=0)

    df_merged.to_csv(out_path+".rfl",header=False,index=False,float_format='%10f', sep=" ", quoting=csv.QUOTE_NONE, escapechar=" ")

    return out_path



if __name__ == "__main__":
    RESULTS_LOCAL_PATH = "results/"
    src_file_path = inspect.getfile(lambda: None)
    scattering_dir = path.abspath(path.join(src_file_path ,"../"))+"/"
    assert len(sys.argv)==3 or len(sys.argv)==4 , "Usage: <scattering_results_handle> <number_of_result_files_to_sample> <(optional) file nametag override>"
    result_handle, num_to_sample = sys.argv[1:3]
    tag_override=None
    if len(sys.argv)==4:
        tag_override = sys.argv[3]
    num_to_sample = int(num_to_sample)
    
   



    reflections_dir =scattering_dir+"random_sample/reflections/"
    sample_results_and_create_scalepack(result_handle, num_to_sample, scattering_dir+RESULTS_LOCAL_PATH, tag_override=tag_override,reflections_dir=reflections_dir,scaling_thing=True)
    
    # reflections_handle = sample_results_and_create_reflection_file(result_handle, num_to_sample, scattering_dir+RESULTS_LOCAL_PATH, tag_override=tag_override,out_directory=reflections_dir)
    # reflections_handle = reflections_handle.split("/")[-1]

    # rfl_to_sca(reflections_handle,reflections_dir=reflections_dir,out_directory=scattering_dir+"random_sample/scalepack/",scaling_method=scaling_method)
    # rfl_to_sca(reflections_handle+"_unmerged",reflections_dir=reflections_dir,out_directory=scattering_dir+"random_sample/scalepack/",scaling_method=scaling_method)