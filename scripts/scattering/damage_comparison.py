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

PDB_PATHS = dict( # contains the target that should be used for the folder.
    neutze =  "/home/speno/AC4DC/scripts/scattering/2lzm.pdb",
    hen = "/home/speno/AC4DC/scripts/scattering/4et8.pdb",
    tetra = "/home/speno/AC4DC/scripts/scattering/5zck.pdb" ,
)            

def multi_damage(batch_handle,params,allowed_atoms_1,CNO_to_N,same_deviations,chosen_folder = None,realistic_crystal_growing_mode = False):
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
    import inspect
    src_file_path = inspect.getfile(lambda: None)
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
            R = stylin(exp_name1,exp_name2,experiment1.max_q,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        else:
            exp1_orientations = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, results_parent_dir=results_batch_dir, **run_params["laser"])
            if exp_name2 != None:
                # sync orientation with damage simulation
                experiment2.set_orientation_set(exp1_orientations)  
                run_params["laser"]["random_orientation"] = False 
                experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, results_parent_dir=results_batch_dir, **run_params["laser"])
            R = stylin(exp_name1,exp_name2,experiment1.max_q,SPI=False,results_parent_dir = results_batch_dir)
        R_data.append([energy[i],fwhm[i],photon_count[i],R]) 
    
    print("R_data:",R_data)
    return np.array(R_data,dtype=np.float64)  

allowed_atoms = ["C_fast","N_fast","O_fast","S_fast"]
CNO_to_N = True
same_deviations = True # whether same position deviations between damaged and undamaged crystal (SPI only)
batch_dir = None # Optional: Specify existing parent folder for batch of results, to add these orientation results to.

R_data = multi_damage("tetra",imaging_params.default_dict,allowed_atoms,CNO_to_N,same_deviations,batch_dir)

##%%
#%%
def plot_that_funky_thing(R_data,cmin=0.1,cmax=0.3):
    df = pd.DataFrame(dict(
        photons =       R_data[:,2],
        fwhm =          R_data[:,1],
        energy =        R_data[:,0], 
        R =             R_data[:,3],
        dummy_column_for_size = R_data[:,0]*0+1
    ))
    print(df)
    # fig = px.scatter_3d(df, x='energy', y='fwhm', z='photons',
    #           color='R', color_continuous_scale=[(0.1, "green"),(0, "green"), (0.2, "yellow"), (0.4, "red"),(1, "red")],range_color=[0,1])   
    # fig = px.scatter_3d(df, x='energy', y='fwhm', z='photons',
    #           color='R',  color_continuous_scale="PuRd", color_continuous_midpoint=0.2) #px.colors.diverging.Armyrose#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
    fig = px.scatter_3d(df, x='photons', y='fwhm', z='energy',
              color='R',  color_continuous_scale="amp", range_color=[cmin,cmax]) #px.colors.diverging.Armyrose#"plasma"#"YlGnBu_r"#cc.m_fire#"inferno"#cmr.ghostlight#cmr.prinsenvlag_r#cmr.eclipse#cc.m_bjy#"viridis"#'Greys'#'binary'
    fig.show()    
    
    unique_E = df.energy.unique()
    #Store data frames
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_E}    
    for key in data_frame_dict.keys():
        data_frame_dict[key] = df[:][df.energy == key]    
    for energy in unique_E:
        df = data_frame_dict[energy] 
        print("Energy:",energy,"eV")
        fig = px.scatter(df,x='photons',y='fwhm',
            color='R',  color_continuous_scale="amp", range_color=[cmin,cmax], size='dummy_column_for_size',size_max=15,)
        fig.show()
        # X,Y,Z = df['photons'].to_numpy(), df['fwhm'].to_numpy(),df['R'].to_numpy()
        # Z_mesh = np.empty((len(X),len(Y)))
        # for i,x in enumerate(X):
        #         Z_mesh[i,i] = Z[i]

        # plt.contourf([X, Y], Z,  alpha=0.7, cmap=plt.cm.jet)
        # fig = go.Figure(data =
        #     go.Contour(
        #         z=[[10, 10.625, 12.5, 15.625, 20],
        #         [5.625, 6.25, 8.125, 11.25, 15.625],
        #         [2.5, 3.125, 5., 8.125, 12.5],
        #         [0.625, 1.25, 3.125, 6.25, 10.625],
        #         [0, 0.625, 2.5, 5.625, 10]],
        #         colorscale='Electric',
        #     ))
        # fig.show()    
plot_that_funky_thing(R_data,0,0.2)
print("Done")

# %%
