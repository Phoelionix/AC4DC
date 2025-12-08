import matplotlib.pyplot as plt
import numpy as np
import sys,os
import matplotlib
import h5py
import os.path as path
from scipy.stats import norm

GUI=True
block=True
separate=False
if not GUI:
    matplotlib.use("Agg")

assert len(sys.argv)>=3
xvg_folders=sys.argv[2:]

centred_on_zero=True

def plot_data(fig:plt.Figure,ax:plt.Axes,x,y,fluence):
    ax.plot(x,y,label=f"{fluence:.1e} J/cm$^{2}$")
    #plt.close()
def plot_gaussian(ax:plt.Axes,x,fwhm, centre):
    sigma = fwhm * np.sqrt(2) / ( np.sqrt(2 * np.log(2)) * 2 )
    gauss = norm.pdf(x,centre,sigma)
    gauss/=max(gauss)
    gauss*=ax.get_ylim()[1]*0.9
    ax.plot(x, gauss,linestyle='dashed',color="dimgrey",label=f"Pulse profile, {fwhm} fs FWHM")

out_tag=sys.argv[1]
duration_fs=10

y_values=[] 
x=None   
fig, ax = plt.subplots(figsize=(6.5,5.5))
fluences=[]
discard_trajectories_with_less_steps=True
for folder in xvg_folders:
    runid = int(path.basename(path.dirname(folder)))
    with h5py.File(path.abspath(path.join(__file__ ,"../../../scripts/molDStructConversion/")) +f"/hdf5_files/{duration_fs}fs.h5", "r") as f:
        fluence = f["fluence"][runid]
    print(folder)
    xvg_file_paths = [path.join(folder,file) for file in os.listdir(folder) if file[-4:]==".xvg"]
    print(sorted(file for file in os.listdir(folder) if file[-4:]==".xvg"))
    if len(xvg_file_paths)==0:
        continue
    for f in xvg_file_paths:
        new_x,y = np.loadtxt(f,comments=("@","#"),unpack=True)
        new_x*=1e3 # ps to fs
        if centred_on_zero:
            new_x = new_x-(max(new_x)-min(new_x))/2
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


plot_gaussian(ax,x,10,(min(x)+max(x))/2)
#ax.set_xlim((max(x)/2,max(x)))
ax.set_xlim((min(x),max(x)))
ax.set_ylim((0,None))
ax.set_ylabel("rmsd $(\\mathrm{\\AA})$")
ax.set_xlabel("time (fs)")

#fig.suptitle()
handles, labels = ax.get_legend_handles_labels()
order = list(reversed(np.argsort([int(fl) for fl in fluences])))
order = order + [i for i in range(len(labels)) if i not in order]
ax.legend(np.array(handles)[order],np.array(labels)[order], loc="upper left")


save_path=f"{path.abspath(path.join(__file__ ,'../../'))}/output/rms_mean_{duration_fs}fs-{out_tag}.png"
fig.savefig(save_path)
if GUI:
    plt.show(block=False)
    plt.pause(0.5)

if block and GUI:
    plt.show()