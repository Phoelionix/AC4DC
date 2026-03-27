#%%
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import matplotlib


folder = "/home/speno/AC4DC/scripts/scattering/SPI_out/giant_hemoglobin_single_traj_data/"
num_bins=10


resolution_range = [60, 2.066] # 2.066 = highest res for 6 keV 
all_rmsds=[]
all_mean_I=[]
all_I_by_time=[]
all_rmsd_by_time=[]
for file in os.listdir(folder):
    data=[]
    with open(os.path.join(folder,file)) as f:
        for line in f.readlines()[1:]:
            entries=line.split(', ')
            intensities =  [float(e) for e in entries[3:]]
            resolution=float(entries[2])
            if not (resolution_range[1] <= resolution <= resolution_range[0]):
                continue

            
            data.append((resolution,intensities))
    
    data.sort(key=lambda resI: resI[0])

    num_per_bin=int(np.ceil(len(data)/num_bins)) 
    #print(len(data))  
    #print([i*num_per_bin for i in range(num_bins)])
    bin_edges = [data[i*num_per_bin][0] for i in range(num_bins)] + [data[-1][0],]  

    
    
    
    binned_intensities = []
    data_resolutions = [entry[0] for entry in data]
    last_idx = 0
    for i in range(len(bin_edges)-1):
        next_idx=data_resolutions.index(bin_edges[i+1])
        binned_intensities.append([v[1] for v in data[last_idx:next_idx]])
        last_idx=next_idx

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

snapshot=None # None 0 6
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

    
print(rmsd_vals.shape)
print(I_vals.shape)
assert rmsd_vals.shape==I_vals.shape, (rmsd_vals.shape, I_vals.shape)

col1=matplotlib.cm.get_cmap('tab10')(0)
col2=matplotlib.cm.get_cmap('tab10')(3)
ax1 = plt.gca()
#ax1.stairs(rmsd_vals,range(len(rmsd_vals)+1),color=col1); ax1.set_ylabel("$I_\sigma$",color=col1) 
ax1.stairs(rmsd_vals/I_vals,range(len(rmsd_vals)+1),color=col1); ax1.set_ylabel("$<I>/I_\sigma$") 
ax1.set_xlabel("Resolution (Å)")
plt.xticks(ticks=range(len(rmsd_vals)+1),labels=[f'{b:.2f}' for b in bin_edges])

ax1.set_xlim(0,len(rmsd_vals))
ax1.invert_xaxis()
#ax.plot(range(len(rmsd_vals)),np.swapaxes(all_mean_I,0,1),)


#ax1.spines['left'].set_color(col2)
#ax1.tick_params(axis='y', colors=col1)
#ax2.spines['right'].set_color(col2)
#ax2.set_ylim((0,ax2.get_ylim()[1]))

if snapshot is None:
    plt.gca().set_title("Integrated")
else:
    snapshot_times = {0:-22.5, 6:22.5}# TEMPORARY
    plt.gca().set_title(f"Snapshot {snapshot_times[snapshot]} fs")


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
