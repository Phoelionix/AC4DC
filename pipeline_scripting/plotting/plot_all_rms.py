import matplotlib.pyplot as plt
import numpy as np
import sys,os
import matplotlib
import h5py
import os.path as path

GUI=True
block=True
separate=False
if not GUI:
    matplotlib.use("Agg")

assert len(sys.argv)>=3
xvg_folders=sys.argv[2:]



def plot_data(fig:plt.Figure,ax:plt.Axes,x,y,fluence):
    ax.plot(x,y,label=f"{fluence:.3e} (J/cm$^{2}$)")
    #plt.close()



out_tag=sys.argv[1]
duration_fs=10

y_values=[] 
x=None   
fig, ax = plt.subplots()
fluences=[]
discard_trajectories_with_less_steps=True
for folder in xvg_folders:
    runid = int(path.basename(path.dirname(folder)))
    with h5py.File(path.abspath(path.join(__file__ ,"../../../scripts/molDStructConversion/")) +f"/hdf5_files/{duration_fs}fs.h5", "r") as f:
        fluence = f["fluence"][runid]
    print(folder)
    print(list(os.listdir(folder)))
    xvg_file_paths = [path.join(folder,file) for file in os.listdir(folder) if file[-4:]==".xvg"]
    if len(xvg_file_paths)==0:
        continue
    for f in xvg_file_paths:
        new_x,y = np.loadtxt(f,comments=("@","#"),unpack=True)
        if x is None:
            x=new_x
        if discard_trajectories_with_less_steps:
            if len(new_x)<len(x):
                continue
            if len(new_x)>len(x):
                y_values=[]
                x=new_x
        assert np.all(x==new_x) 
        y_values.append(y)    

    y_mean = np.mean(y_values,axis=0)
    plot_data(fig,ax,x,y_mean, fluence=fluence)
    fluences.append(fluence)

#fig.suptitle()
handles, labels = plt.gca().get_legend_handles_labels()
order = list(reversed(np.argsort([int(fl) for fl in fluences])))
fig.legend(np.array(handles)[order],np.array(labels)[order])


save_path=f"{path.abspath(path.join(__file__ ,'../../'))}/output/rms_mean-{out_tag}.png"
fig.savefig(save_path)
if GUI:
    plt.show(block=False)
    plt.pause(0.5)

if block and GUI:
    plt.show()