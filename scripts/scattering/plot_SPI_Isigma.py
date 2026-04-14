#%%
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
sys.path.append('/home/speno/AC4DC/scripts/')
from plotter_core import Plotter
from scatter import res_to_q, q_to_res

SNAPSHOTS_NEED_PULSE_I_SCALING=True; plasma_sim_handle="giant_hemoglobin_solvated_1"
EVEN_BIN_ALLOCATION=False
#folder = "/home/speno/AC4DC/scripts/scattering/SPI_out/giant_hemoglobin_single_traj_data/"
folder = "/home/speno/AC4DC/scripts/scattering/SPI_out/giant_hemoglobin_single_pattern_500x500px/"
num_bins=50


TEMP_TIMES = np.array([0.0075, 0.015, 0.0225, 0.03, 0.0375, 0.045, 0.0525])*1e3-30   


if SNAPSHOTS_NEED_PULSE_I_SCALING:
    pl = Plotter(plasma_sim_handle)



resolution_range = [9999, 2.066] # 2.066 = highest res for 6 keV 
all_rmsds=[]
all_mean_I=[]
all_I_by_time=[]
all_rmsd_by_time=[]
for file in list(os.listdir(folder)[:2]):
    data=[]
    with open(os.path.join(folder,file)) as f:

        if SNAPSHOTS_NEED_PULSE_I_SCALING:
            times=TEMP_TIMES # FIXME read this from header
            t_idx, t =pl.get_nearest_time(TEMP_TIMES)
            pulse_intensity_arr = pl.intensityData[t_idx]
        for line in f.readlines()[1:]:
            entries=line.split(', ')
            intensities =  np.array([float(e) for e in entries[3:]])
            if SNAPSHOTS_NEED_PULSE_I_SCALING:
                intensities*=pulse_intensity_arr
            resolution=float(entries[2])
            if not (resolution_range[1] <= resolution <= resolution_range[0]):
                continue

            
            data.append((resolution,intensities))
    
    data.sort(key=lambda resI: resI[0])

    data_resolutions = [entry[0] for entry in data]

    if EVEN_BIN_ALLOCATION:
        num_per_bin=int(np.ceil(len(data)/num_bins)) 
        print(f"Aiming for {num_per_bin} reflections per bin")
        #print(len(data))  
        #print([i*num_per_bin for i in range(num_bins)])
        bin_edges = [data[i*num_per_bin][0] for i in range(num_bins)] + [data[-1][0],]  
    else: # even spacing of q
        bin_edges=np.linspace(res_to_q(data[0][0]),res_to_q(data[-1][0]),num_bins+1,endpoint=True)
        bin_edges = q_to_res(np.array(bin_edges))
        bin_edges = [data_resolutions[idx] for idx in np.searchsorted(data_resolutions,bin_edges) ]
        #bin_edges = list(sorted(bin_edges))

    binned_intensities = []
    last_idx = 0
    for i in range(len(bin_edges)-1):
        next_idx=data_resolutions.index(bin_edges[i+1])
        binned_intensities.append([v[1] for v in data[last_idx:next_idx]])
        last_idx=next_idx
        if EVEN_BIN_ALLOCATION:
            if abs(len(binned_intensities[-1])-num_per_bin) > 0.2*num_per_bin:
                print(f"Warning, bin {i}, {bin_edges[i]} - {bin_edges[i+1]} A has: {len(binned_intensities)[-1]}")

    rmsd_array=[]
    mean_I_array=[]
    all_rmsds.append(rmsd_array)
    all_mean_I.append(mean_I_array)
    
    r_arr=[]; I_arr=[]
    all_rmsd_by_time.append(r_arr)
    all_I_by_time.append(I_arr)
    for i, intensities in enumerate(binned_intensities):
        bin = bin_edges[i:i+2]
        tot_I = np.sum(np.array(list(intensities),dtype=float),axis=-1)
        mean_I = np.mean(np.array(tot_I))
        rmsd = np.sqrt(np.mean([np.sum((I-mean_I)**2) for I in tot_I]))
        rmsd_array.append(rmsd)
        mean_I_array.append(mean_I)

        mini_rmsd=[]
        mini_I=[]
        r_arr.append(mini_rmsd)
        I_arr.append(mini_I)
        for I_snapshot in np.array(np.swapaxes(intensities,0,1)):
            mini_rmsd.append(np.sqrt(np.mean([np.sum((I-np.mean(I_snapshot))**2) for I in I_snapshot])))
            mini_I.append(np.mean(I_snapshot))
        # print(f"bin {bin}; RMSD {rmsd:.2e}, mean intensity: {mean_I:.2e} mean intensity each snapshot:  "
        #       + ', '.join([f"{v:.2e}" for v in np.mean(intensities,axis=0)]))
    
#%%

stairs_plot=True
fig,ax1 = plt.subplots(figsize=(9,4))
ax2 = ax1.twiny()
#for snapshot in [None]:
#for snapshot in [0,1,2,3,4,5,6]:
#for snapshot in [0,6]:
for snapshot, label in zip([None,0,6],["Integrated","$t \\;\\; = -22.5$ fs","$t \\;\\; = +22.5$ fs"]):
    if snapshot is None:
        _rmsd_arr=all_rmsds
        _I_arr=all_mean_I
    else:
        print(np.array(all_rmsd_by_time,dtype=object).shape)
        print(np.array(all_I_by_time,dtype=object).shape)
        _rmsd_arr=np.array(all_rmsd_by_time,dtype=float)[...,snapshot]
        _I_arr=np.array(all_I_by_time,dtype=float)[...,snapshot]
        print(_rmsd_arr.shape)
        print(_I_arr.shape)
    rmsd_vals=np.mean(np.array(_rmsd_arr),axis=0)
    I_vals=np.mean(np.array(_I_arr),axis=0)
    bin_edges=np.array(bin_edges)
        
    print(rmsd_vals.shape)
    print(I_vals.shape)
    assert rmsd_vals.shape==I_vals.shape, (rmsd_vals.shape, I_vals.shape)

    col1=matplotlib.cm.get_cmap('tab10')(0)
    col2=matplotlib.cm.get_cmap('tab10')(3)
    #ax1 = plt.gca()
    #ax1.stairs(rmsd_vals,range(len(rmsd_vals)+1),color=col1); ax1.set_ylabel("$I_\sigma$",color=col1) 
    Y=rmsd_vals/I_vals
    if stairs_plot:
        ax1.stairs(Y,range(len(rmsd_vals)+1),label=label,linewidth=2)
        ax1.set_xlim(len(rmsd_vals),0)
        ax2.set_xlim(len(rmsd_vals),0)
        interval=3 # XXX
        ax1.set_xticks(ticks=list(range(len(rmsd_vals)+1))[::interval],
                       #labels=[f'{res_to_q(b):.2f}' for b in bin_edges][::10])
                       labels=[f'{res_to_q(b):.2f}' for b in bin_edges][::interval])
        ax1.set_ylim(0.4,2)
        ax2.set_xticks(ticks=list(range(len(rmsd_vals)+1))[::interval], 
                       labels=[f'{b:.2f}' for b in bin_edges][::interval])

        ax1.hlines(1,colors="black",alpha=1,lw=1,linestyles=[(0, (1, 3))],xmin=len(rmsd_vals),xmax=0)

    #ax.plot(range(len(rmsd_vals)),np.swapaxes(all_mean_I,0,1),)
    else: # points
        ax1.plot((bin_edges[:-1]+bin_edges[1:])/2,Y,marker='o',label=label)
        ax1.set_xlim(bin_edges[-2],0)
        ax1.set_ylim(None,max(Y[:-1]*1.02))
    ax1.set_ylabel("$<I>/\sigma(I)$") 
    ax1.set_xlabel("$q$ (Å$^{-1}$)")
    ax2.set_xlabel("Resolution (Å)")
    ax1.tick_params("x",rotation=90)
    ax2.tick_params("x",rotation=90)
    ax1.legend()

    #ax1.spines['left'].set_color(col2)
    #ax1.tick_params(axis='y', colors=col1)
    #ax2.spines['right'].set_color(col2)
    #ax2.set_ylim((0,ax2.get_ylim()[1]))

    # if snapshot is None:
    #     plt.gca().set_title("Integrated")
    # else:
    #     snapshot_times = {0:-22.5, 6:22.5}# TEMPORARY
    #     plt.gca().set_title(f"Snapshot {TEMP_TIMES[snapshot]} fs")
plt.show()

# print(rmsd_vals.shape)
# print(I_vals.shape)
# assert rmsd_vals.shape==I_vals.shape, (rmsd_vals.shape, I_vals.shape)

# col1=matplotlib.cm.get_cmap('tab10')(0)
# col2=matplotlib.cm.get_cmap('tab10')(3)
# ax1 = plt.gca()
# plt.stairs(rmsd_vals,range(len(rmsd_vals)+1),color=col1)
# ax1.set_xlabel("Resolution (Å)")
# plt.xticks(ticks=range(len(rmsd_vals)+1),labels=[f'{b:.2f}' for b in bin_edges])

# ax1.set_xlim(0,len(rmsd_vals))
# ax2 = ax1.twinx()
# #ax.plot(range(len(rmsd_vals)),np.swapaxes(all_mean_I,0,1),)

# ax2.plot(np.array(range(len(rmsd_vals)))+0.5,I_vals,color=col2)
# ax1.set_ylabel("$I_\sigma$",color=col1) 
# #ax1.spines['left'].set_color(col2)
# ax1.tick_params(axis='y', colors=col1)
# ax2.set_ylabel("$<I>)$",color=col2)
# #ax2.spines['right'].set_color(col2)
# ax2.tick_params(axis='y', colors=col2)
# #ax2.set_ylim((0,ax2.get_ylim()[1]))

# if snapshot is None:
#     plt.gca().set_title("Integrated")
# else:
#     snapshot_times = {0:-22.5, 6:22.5}# TEMPORARY
#     plt.gca().set_title(f"Snapshot {snapshot_times[snapshot]} fs")

            


            







# %%
