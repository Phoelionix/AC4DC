#%%
from scatter import XFEL,Crystal,stylin
import numpy as np
import os.path as path
import os
from shutil import rmtree
import inspect
import copy
import imaging_params
import time
import plotly.graph_objects as go
import pandas as pd
import plotly.express as px
import colorcet as cc; import cmasher as cmr
from core_functions import get_sim_params
import matplotlib.pyplot as plt

import inspect
src_file_path = inspect.getfile(lambda: None)  
my_dir = path.abspath(path.join(src_file_path ,"../")) + "/"

PDB_PATHS = dict( # contains the target that should be used for the folder.
    neutze =  my_dir + "2lzm.pdb",
    hen = my_dir + "4et8.pdb",
    tetra = my_dir+ "5zck.pdb" ,
)          


def multi_damage(batch_handle,params,allowed_atoms_1,CNO_to_N,same_deviations,chosen_folder = None,realistic_crystal_growing_mode = False,get_R_only=True):
    '''
    
    '''
    # Turn off stdv for crystal 
    if params["laser"]["SPI"] == False:
        params["crystal"]["positional_stdv"] = 0
    pdb_path = PDB_PATHS[batch_handle]

    """
    Folder structure:
    ...'root_results_dir'/'results_batch_dir'/<experiment_results>/<orientation>
    """
    sim_data_batch_dir = path.abspath(path.join(src_file_path ,"../../../output/__Molecular/"+batch_handle)) + "/"
    sim_input_dir = path.abspath(path.join(src_file_path ,"../../../input/")) + "/"
    # Folder containing the folders corresponding to each batch of handles.
    root_results_dir = path.join(src_file_path,"results_R_plotting/")+ "/"
    if path.isdir(root_results_dir):
        rmtree(root_results_dir)
    if path.isfile(root_results_dir):
        raise Exception("target directory is existing file's name.")
    
    # Folder containing the batch of handles (i.e. the results' folders' parent directory)
    results_batch_dir = chosen_folder
    if chosen_folder is None:
        version_number = 1
        count = 0
        while True:
            if count > 99:
                raise Exception("could not find valid parent folder name in " + str(count) + " loops")
            results_batch_dir = root_results_dir+batch_handle +"_v"+str(version_number) +"/" # needs to be synced with other functions
            if path.exists(path.dirname(results_batch_dir)):
                version_number+=1
                count+=1
                continue 
            break
    os.makedirs(results_batch_dir)
    # Run imagings for each simulation 
    exp1_qualifier = "_real"
    exp2_qualifier = "_ideal"        

    # find shape of array of params
    num_results = len(os.listdir(sim_data_batch_dir))
    sim_handle_list, start_time, end_time, energy, fwhm, photon_count = ([None]*num_results for i in range(6))
    for i, sim_handle in enumerate(os.listdir(sim_data_batch_dir)):
        sim_handle_list[i] = sim_handle
        start_time[i],end_time[i],energy[i],fwhm[i],photon_count[i],param_dict,unit_dict = get_sim_params(sim_input_dir,sim_data_batch_dir,sim_handle)
    
    R_data = []
    for i, sim_handle in enumerate(sim_handle_list):
        run_params = copy.deepcopy(params)
        
        #TODO print properly
        print(sim_handle)
        print(energy[i],fwhm[i],photon_count[i])
        print(param_dict)
        print(unit_dict)
        exp_name1 = sim_handle + "_" + exp1_qualifier
        exp_name2 = sim_handle + "_" + exp2_qualifier


        # Spend a few billion on XFELs
        experiment1 = XFEL(exp_name1,energy[i],**run_params["experiment"])
        experiment2 = XFEL(exp_name2,energy[i],**run_params["experiment"])
        # Wait a few months for crystals to grow

        # For added realism
        if realistic_crystal_growing_mode:
            for day in range(90):
                print("Crystals growing,",90-day,"days remaining...")
                time.sleep(216000)
        else:
            print("Crystals growing at breakneck speeds...")
            time.sleep(0.3)            

        crystal = Crystal(pdb_path,allowed_atoms_1,is_damaged=True,CNO_to_N = CNO_to_N, **run_params["crystal"])
        # The undamaged crystal uses the initial state but still performs the same integration step with the pulse profile weighting.
        # we copy the other crystal so that it has the same deviations in coords
        if same_deviations:
            # we copy the other crystal so that it has the same deviations in coords
            crystal_undmged = copy.deepcopy(crystal)
            crystal_undmged.is_damaged = False
        else:
            crystal_undmged = Crystal(pdb_path,allowed_atoms_1,is_damaged=False,CNO_to_N = CNO_to_N, **run_params["crystal"])
        if i == 0:
            crystal.plot_me(250000)
    
        if params["laser"]["SPI"]:
            SPI_result1 = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal,results_parent_dir=results_batch_dir, **run_params["laser"])
            SPI_result2 = experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged,results_parent_dir=results_batch_dir, **run_params["laser"])
            R = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        else:
            exp1_orientations = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, results_parent_dir=results_batch_dir, **run_params["laser"])
            if exp_name2 != None:
                # sync orientation with damage simulation
                experiment2.set_orientation_set(exp1_orientations)  
                run_params["laser"]["random_orientation"] = False 
                experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, results_parent_dir=results_batch_dir, **run_params["laser"])
            R = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=False,results_parent_dir = results_batch_dir)
        R_data.append([energy[i],fwhm[i],photon_count[i],R]) 
    
    print("R_data:",R_data)
    return np.array(R_data,dtype=np.float64)  

allowed_atoms = ["C_fast","N_fast","O_fast","S_fast"]
CNO_to_N = True
same_deviations = True # whether same position deviations between damaged and undamaged crystal (SPI only)
batch_dir = None # Optional: Specify existing parent folder for batch of results, to add these orientation results to.

R_data_crys = multi_damage("tetra",imaging_params.tetra_dict,allowed_atoms,CNO_to_N,same_deviations,batch_dir,get_R_only=True)
R_data_SPI = multi_damage("tetra",imaging_params.tetra_dict_SPI,allowed_atoms,CNO_to_N,same_deviations,batch_dir,get_R_only=True)

#%%
def plot_that_funky_thing(R_data,cmin=0.1,cmax=0.3,clr_scale="amp",**kwargs):
    log_photons = np.log(R_data[:,2])-np.min(np.log(R_data[:,2]))
    log_photons += np.max(log_photons)/2
    # photon_size_thing = R_data[:,2] + np.max(R_data[:,2])/10
    photon_size_thing = 1+np.square(np.log(R_data[:,2])-np.min(np.log((R_data[:,2]))))

    df = pd.DataFrame(dict(
        photons =       R_data[:,2],
        log_photons = log_photons,
        photon_size_thing = photon_size_thing,
        fwhm =          R_data[:,1],
        energy =        R_data[:,0], 
        R =             R_data[:,3],
        dummy_column_for_size = R_data[:,0]*0+1
    ))
    print(df)
    # fig = px.scatter_3d(df, x='energy', y='fwhm', z='photons',
    #           color='R', color_continuous_scale=[(0.1, "green"),(0, "green"), (0.2, "yellow"), (0.4, "red"),(1, "red")],range_color=[0,1])   
    # fig = px.scatter_3d(df, x='energy', y='fwhm', z='photons',
    #           color='R',  color_continuous_scale="PuRd", color_continuous_midpoint=0.2) 
    # fig = px.scatter_3d(df, x='photons', y='fwhm', z='energy',
    #           color='R',  color_continuous_scale="amp", range_color=[cmin,cmax])'
    fig = px.scatter_3d(df, x='photons', y='fwhm', z='energy',
              color='R',  color_continuous_scale=clr_scale, range_color=[cmin,cmax],**kwargs)
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.25, y=-1.25, z=1.25)
    )

    fig.update_layout(scene_camera=camera)         
    fig.show()    

    fig = px.scatter(df, x='fwhm', y='R', symbol='energy', size='log_photons',
              color='R', color_continuous_scale=clr_scale, range_color=[cmin,cmax],size_max=15,opacity=1,**kwargs) 
    fig.update_layout(
        showlegend=True,
        legend=dict(
        yanchor="top",
        y=0.98,
        xanchor="left",
        x=1.13,
    ))      
    fig.show()

    # Dot size represents num photons, 
    df.sort_values('photon_size_thing')
    fig = px.scatter(df, x='fwhm', y='energy', size='photon_size_thing',
              color='R', color_continuous_scale=clr_scale, range_color=[cmin,cmax],size_max=30,opacity=1,**kwargs)     
    fig.show()    
        
    # fig = px.scatter(df, x='fwhm', y='R', symbol='fwhm', size='log_photons',
    #           color='energy', color_continuous_scale = 'plasma',size_max=15,opacity=1,**kwargs) 
    # fig.update_layout(
    #     showlegend=True,
    #     legend=dict(
    #     yanchor="top",
    #     y=0.98,
    #     xanchor="left",
    #     x=1.13,
    # ))     
    # fig.show()
    # fig = px.scatter(df, x='fwhm', y='R', size='log_photons',
    #           color='energy', color_continuous_scale = 'plasma',size_max=15,opacity=1,**kwargs) 
    # fig.show()    



    ''' plot separate energies
    # Get min and max for each axis for consistency
    ranges = {}
    for col in df:
        range = df[col].agg(['min','max'])
        range[0]-=range[1]/10 
        range[1]+= range[1]/10
        ranges[col] = range

    #get data frames for each energy
    unique_E = np.sort(df.energy.unique())
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_E}    
    for key in data_frame_dict.keys():
        data_frame_dict[key] = df[:][df.energy == key]    
    for energy in unique_E:
        df = data_frame_dict[energy] 
        print("Energy:",energy,"eV")
        fig = px.scatter(df,x='photons',y='fwhm',
            color='R',  color_continuous_scale=clr_scale, range_color=[cmin,cmax], size='dummy_column_for_size',size_max=15,opacity=1,**kwargs)
        fig.update_xaxes(range=ranges['photons'])
        fig.update_yaxes(range=ranges['fwhm'])           
        fig.show()
    '''

print("-----------------Crystal----------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
plot_that_funky_thing(R_data_crys,0.05,0.25,"temps",template="plotly_dark") #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
print("-------------------SPI------------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
plot_that_funky_thing(R_data_SPI,0.05,0.25,"temps",template="plotly_dark") #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
print("Done")

#TODO store results.

# %%
