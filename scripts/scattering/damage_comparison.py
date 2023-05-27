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

# full sim that does point sampling, thus assumes a "smooth" distribution. Thus good for low cell count input. Num cells is entire whole target.
# Compared with non-SPI (crystal) which assumes crystal and does bragg points. Num cells is supercell.
#TODO rename the diff sims, 'non-SPI' unclear.
#TODO currently all cells are used for calculation of R, but we may wish to ignore cells outside the considered radius? Though I haven't see anyone do this, it would depend on what image reconstruction algos do I guess. But we get empty corners with high res. limits so


MODE_DICT = {0:'crystal',1:'spi',2:'both'}
NO_CRYSTAL_DEV = False

PDB_PATHS = dict( # <value:> the target that should be used for <key:> the name of simulation output batch folder.
    tetra = my_dir+ "targets/5zck.pdb",
    lys = my_dir + "targets/4et8.pdb", #"targets/2lzm.pdb",
    test = my_dir + "targets/4et8.pdb", 
    lys_tmp = my_dir + "targets/4et8.pdb",
    
    fcc = "targets/FCC.pdb",
)          


def multi_damage(params,pdb_path,allowed_atoms_1,CNO_to_N,S_to_N,same_deviations,plasma_batch_handle = "", plasma_handles = None, sctr_results_batch_dir = None, get_R_only=True, realistic_crystal_growing_mode = False,specific_energy = None):
    '''
    NB: "Plasma" is used a shortened way to refer to the output of the damage simulation  
    parameters:

    If plasma_batch_handle is not provided, handle_list must be provided to prevent.
    plasma_batch_handle<optional>:  name of folder containing individual simulation output folders. If not provided, with use handle_list in the default directory.
    handle_list <optional>: list of folder names to iterate through within the batch folder. If not provided, all are used. 
    '''
    if plasma_batch_handle == "":
        assert plasma_handles is not None

    # Turn off stdv for crystal  #TODO tackle this better
    if NO_CRYSTAL_DEV and params["laser"]["SPI"] == False:
        params["crystal"]["positional_stdv"] = 0

    if params["laser"]["SPI"] == False:  
        # Folder structure:
        #.../'root_results_dir'/'sctr_results_batch_dir'/<experiment_results>/<orientation>
        # Folder containing the folders corresponding to each batch of handles.
        root_results_dir = path.join(my_dir,"results_R_plotting/")+ "/"
        if path.isdir(root_results_dir):
            rmtree(root_results_dir)
        if path.isfile(root_results_dir):
            raise Exception("target directory is existing file's name.")
       
    
        # Determine the scatter results' folders' parent directory. (used for infinite crystal only)
        if sctr_results_batch_dir is None:
            version_number = 1
            count = 0
            while True:
                if count > 99:
                    raise Exception("could not find valid parent folder name in " + str(count) + " loops")
                sctr_results_batch_dir = root_results_dir+plasma_batch_handle +"_v"+str(version_number) +"/" # needs to be synced with other functions
                if path.exists(path.dirname(sctr_results_batch_dir)):
                    version_number+=1
                    count+=1
                    continue 
                break
        os.makedirs(sctr_results_batch_dir)
    # Run imagings for each simulation 
    exp1_qualifier = "_real"
    exp2_qualifier = "_ideal"        

    sim_data_batch_dir = path.abspath(path.join(my_dir ,"../../output/__Molecular/"+plasma_batch_handle)) + "/"
    sim_input_dir = path.abspath(path.join(my_dir ,"../../input/")) + "/"
    if plasma_handles == None:
        plasma_handles = []
        print("Iterating through entire batch of plasma simulation results.")
        for i, sim_handle in enumerate(os.listdir(sim_data_batch_dir)):
            plasma_handles.append(sim_handle)
    else:
        print("Iterating through",plasma_handles)

    start_time, end_time, energy, fwhm, photon_count = ([None]*len(plasma_handles) for _ in range(5))
    necessary_files = ["freeDist.csv","intensity.csv"]
    for i, sim_handle in enumerate(plasma_handles):
        # Check some of the expected files are present (not checking for existence of bound state files at present.) 
        plasma_output = sim_data_batch_dir + sim_handle
        assert path.isdir(plasma_output), "Directory not found for "+plasma_output
        for dat_file in necessary_files: 
            assert path.isfile(plasma_output + "/" +dat_file), "Missing "+dat_file+" in "+plasma_output
        # Index parameters of simulation
        start_time[i],end_time[i],energy[i],fwhm[i],photon_count[i],param_dict,unit_dict = get_sim_params(sim_input_dir,sim_data_batch_dir,sim_handle)
    R_data = []
    names = []
    for i, sim_handle in enumerate(plasma_handles):
        if specific_energy is not None and energy[i] != specific_energy:
            continue
        names.append(sim_handle)
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
        if realistic_crystal_growing_mode:
            # For added realism
            for day in range(90):
                print("Crystals growing,",90-day,"days remaining...")
                time.sleep(216000)
        else:
            print("Crystals growing at breakneck speeds...")
            time.sleep(0.3)            

        # The undamaged crystal form factor is constant (the initial state's) but still performs the same integration step with the pulse profile weighting.
        crystal = Crystal(pdb_path,allowed_atoms_1,is_damaged=True,CNO_to_N = CNO_to_N, S_to_N = S_to_N, **run_params["crystal"])
        if same_deviations:
            # we copy the other crystal so that it has the same deviations in coords
            crystal_undmged = copy.deepcopy(crystal)
            crystal_undmged.is_damaged = False
        else:
            crystal_undmged = Crystal(pdb_path,allowed_atoms_1,is_damaged=False,CNO_to_N = CNO_to_N,S_to_N = S_to_N, **run_params["crystal"])
        if i == 0:
            crystal.plot_me(250000)
    
        if params["laser"]["SPI"]:
            SPI_result1 = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, **run_params["laser"])
            SPI_result2 = experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, **run_params["laser"])
            R = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=True,SPI_max_q = None,SPI_result1=SPI_result1,SPI_result2=SPI_result2)
        else:
            exp1_orientations = experiment1.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            if exp_name2 != None:
                # sync orientation with damage simulation
                experiment2.set_orientation_set(exp1_orientations)  
                run_params["laser"]["random_orientation"] = False 
                experiment2.spooky_laser(start_time[i],end_time[i],sim_handle,sim_data_batch_dir,crystal_undmged, results_parent_dir=sctr_results_batch_dir, **run_params["laser"])
            R = stylin(exp_name1,exp_name2,experiment1.max_q,get_R_only=get_R_only,SPI=False,results_parent_dir = sctr_results_batch_dir)
        R_data.append([energy[i],fwhm[i],photon_count[i],R]) 
    
    print("R_data:",R_data)
    return np.array(R_data,dtype=np.float64), names, param_dict

def plot_that_funky_thing(R_data,names,param_dict,cmin=0.1,cmax=0.3,clr_scale="amp",**kwargs):
    log_photons = np.log(R_data[:,2])-np.min(np.log(R_data[:,2]))
    log_photons += np.max(log_photons)/2
    # photon_size_thing = R_data[:,2] + np.max(R_data[:,2])/10
    photon_size_thing = 1+np.square(np.log(R_data[:,2])-np.min(np.log((R_data[:,2]))))

    assert param_dict[3] == "R"
    photon_measure = param_dict[2]
    df = pd.DataFrame({
        "name": names, 
        "energy": R_data[:,0],
        "fwhm": R_data[:,1],
        photon_measure: R_data[:,2],
        "R": R_data[:,3],
        #"log_photons": log_photons,
        "photon_size_thing": photon_size_thing, 
        "_": R_data[:,0]*0+1, # dummy column for size.
    })
    print(df)
    #3D
    fig = px.scatter_3d(df, x=photon_measure, y='fwhm', z='energy',
              color='R',  color_continuous_scale=clr_scale, range_color=[cmin,cmax],**kwargs)
    camera = dict(
        up=dict(x=0, y=0, z=1),
        center=dict(x=0, y=0, z=0),
        eye=dict(x=1.25, y=-1.25, z=1.25)
    )

    fig.update_layout(scene_camera=camera)         
    fig.show()    

    label_col = 'name' 
    if len(df['name']) > 10:
        print("Turning off labels due to high number of points")
        label_col = None    
    fig = px.scatter(df, text=label_col,x='fwhm', y='R', symbol='energy', size='photon_size_thing',#'log_photons',
              color='R', color_continuous_scale=clr_scale, range_color=[cmin,cmax],size_max=15,opacity=1,**kwargs) 
    fig.update_layout(
        showlegend=True,
        legend=dict(
        yanchor="top",
        y=0.98,
        xanchor="left",
        x=1.13,
    ))      
    fig.update_traces(textposition='top center')
    fig.show()

    # 'intuitive' plot, where:
    # Dot size represents num photons
    # colour is R (damage)
    # horiz axis is fwhm
    # vert axis is energy of photons 
    df = df.sort_values('photon_size_thing',ascending=False)
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



    #plot separate energies
    # Get min and max to use for axes here, so we can be consistent.
    ranges = {}
    for col in df:
        if col == "name":
            continue
        #rng = df[col].agg(['min','max'])
        rng = [df[col].min(),df[col].max()]
        rng[0]-=rng[1]/10 
        rng[0] = max(0,rng[0])
        rng[1]+= rng[1]/10
        ranges[col] = rng

    original_df = df
    #get data frames for each energy
    unique_E = np.sort(df.energy.unique())
    unique_fwhm = np.sort(df.fwhm.unique())
    unique_photons = np.sort(df[photon_measure].unique())
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_E}    
    for key in data_frame_dict.keys():
        data_frame_dict[key] = df[:][df.energy == key]    
    for energy in unique_E:
        df = data_frame_dict[energy] 
        print("Energy:",energy,"eV")
        fig = px.scatter(df,x='fwhm',y=photon_measure,
            color='R',  color_continuous_scale=clr_scale, range_color=[cmin,cmax], size='_',size_max=15,opacity=1,**kwargs)
        fig.update_xaxes(range=ranges['fwhm'])
        fig.update_yaxes(type="log",range=np.log10(ranges[photon_measure]))
        fig.show()
    # and do a contour for 12keV
    if 12000 in data_frame_dict.keys(): 
        df = data_frame_dict[12000]
        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df["R"],
                x = df["fwhm"],
                y = [1e12*p for p in df[photon_measure]],

                colorscale = 'amp',#clr_scale, #'electric',
                line_smoothing=0,
                connectgaps = True,
                zmin = 0, zmax = 0.4,
                contours = go.contour.Contours(start = 0.05, end= ranges['R'][1], size = 0.05),
                colorbar = dict(title = "R", tickvals = np.arange(0.05,ranges['R'][1],0.05),)
            ))
        fig.update_xaxes(title="FWHM (fs)",type="log")
        fig.update_yaxes(title="Photons")
        fig.update_yaxes(type="log",range=np.log10(ranges[photon_measure]),tickvals = [1e10,1e11,1e12,1e13],tickformat = '.0e')
        fig.update_layout(width = 750, height = 600,)
        fig.show()
        df = original_df

    #get data frames for 1e12 photon count and do a contour
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_photons}   
    if 1 in data_frame_dict.keys():  
        for key in data_frame_dict.keys():
            data_frame_dict[key] = df[:][df.photons == key]    
        df = data_frame_dict[1]

        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df["R"],
                x = df["fwhm"],
                y = df["energy"],

                colorscale = 'amp',#clr_scale, #'electric',
                line_smoothing=0,
                connectgaps = True,
                zmin = 0, zmax = 0.4,
                contours = go.contour.Contours(start = 0.05, end= ranges['R'][1], size = 0.05),
                colorbar = dict(title = "R", tickvals = np.arange(0.05,ranges['R'][1],0.05),)
            ))
        fig.update_layout(width = 750, height = 600,)
        fig.update_xaxes(title="FWHM (fs)",type="log")
        fig.update_yaxes(title="Energy (eV)")

        fig.show()   
        df = original_df

    #get data frames for 25 fs fwhm and do a contour
    data_frame_dict = {elem : pd.DataFrame() for elem in unique_fwhm}    
    if 25 in data_frame_dict.keys(): 
        for key in data_frame_dict.keys():
            data_frame_dict[key] = df[:][df.fwhm == key]    
        df = data_frame_dict[25]

        fig = go.Figure(data =  
            go.Contour(
                #z=R_mesh,
                #x=unique_fwhm, 
                #y=unique_photons,
                z = df["R"],
                x = df["energy"],
                y = [1e12*p for p in df[photon_measure]],

                colorscale = 'amp',#clr_scale, #'electric',
                line_smoothing=0,
                connectgaps = True,
                zmin = 0, zmax = 0.4,
                contours = go.contour.Contours(start = 0.05, end= ranges['R'][1], size = 0.05),
                colorbar = dict(title = "R", tickvals = np.arange(0.05,ranges['R'][1],0.05),)
            ))
        fig.update_layout(width = 750, height = 600,)
        fig.update_xaxes(title="Energy (eV)")
        fig.update_yaxes(title="Photons",type="log",range=np.log10(ranges[photon_measure]),tickvals = [1e10,1e11,1e12,1e13],tickformat = '.0e')
        fig.show()        
        df = original_df
#%%
#------------Get R--------------------

if __name__ == "__main__":
    kwargs = dict(
        plasma_batch_handle = "", 
        plasma_handles = None, 
        sctr_results_batch_dir = None, 
        get_R_only=True,
    )


    batch_mode = False # Just doing this as I want to quickly switch between batches and specific runs.

    mode = 1  #0 -> infinite crystal, 1 -> finite crystal/SPI, 2-> both
    allowed_atoms = ["C_fast","N_fast","O_fast","S_fast"] 
    
    same_deviations = False # whether same position deviations between damaged and undamaged crystal (SPI only)
    
    if batch_mode:
        CNO_to_N = True
        S_to_N = False
        batch_handle = "tetra" 
        kwargs["plasma_batch_handle"] = batch_handle
        batch_dir = None # Optional: Specify existing parent folder for batch of results, to add these orientation results to.
        pdb_path = PDB_PATHS[batch_handle]
    else:
        CNO_to_N = True
        S_to_N = False
        kwargs["plasma_batch_handle"] = ""
        #kwargs["plasma_handles"] = ["lys_nass_3","lys_nass_no_S_1"]    #Note that S_to_N must be true to compare effect on nitrogen R factor. Comparing S_to_N true with nitrogen only sim, then S_to_N false with nitrogen+sulfur sim, let's us compare the true effect of sulfur on R facator.
        kwargs["plasma_handles"] = ["lys_nass_3"]  
        #pdb_path = PDB_PATHS["fcc"]
        pdb_path = PDB_PATHS["lys"]

    if MODE_DICT[mode] != "spi":
        scatter_data = multi_damage(imaging_params.default_dict,pdb_path,allowed_atoms,CNO_to_N,S_to_N,same_deviations,**kwargs)
    if MODE_DICT[mode] != "crystal":
        scatter_data = multi_damage(imaging_params.default_dict_SPI,pdb_path,allowed_atoms,CNO_to_N,S_to_N,same_deviations,**kwargs)


#--------------------------------------
#%%
#------------Plot----------------------
if __name__ == "__main__":
    if MODE_DICT[mode] != "spi":
        print("-----------------Crystal----------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(*scatter_data,0,0.20,"temps",template="plotly_dark") # 'electric' #"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    if MODE_DICT[mode] != "crystal":
        print("-------------------SPI------------------------")                                                    #"plotly" #"simple_white" #"plotly_white" #"plotly_dark"
        plot_that_funky_thing(*scatter_data,0,0.20,"temps",template="plotly_dark") # 'electric'#"fall" #"Temps" #"oxy" #RdYlGn_r #PuOr #PiYg_r #PrGn_r
    print("Done")

#--------------------------------------
#TODO store results.
#TODO use form factor file produced by AC4DC rather than recalc so much.
# %%
