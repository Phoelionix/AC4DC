import matplotlib.pyplot as plt
import numpy as np
import sys,os
import matplotlib

GUI=True
block=True
separate=False
if not GUI:
    matplotlib.use("Agg")

assert len(sys.argv)>=2
xvg_file_paths=sys.argv[1:]



def plot_data(x,y,save_path):
    plt.figure()
    plt.plot(x,y)
    if GUI:
        plt.show(block=False)
        plt.pause(0.5)
    plt.savefig(save_path)
    #plt.close()

if separate:
    for f in xvg_file_paths:
        assert f[-4:]==".xvg"
        x,y = np.loadtxt(f,comments=("@","#"),unpack=True)
        plot_data(x,y,save_path=f"{f[:-4]}_rms.png")
else:
    y_values=[] 
    x=None   
    for f in xvg_file_paths:
        new_x,y = np.loadtxt(f,comments=("@","#"),unpack=True)
        if x is None:
            x=new_x
        assert np.all(x==new_x) 
        y_values.append(y)
    y_mean = np.mean(y_values,axis=0)
    plot_data(x,y_mean,save_path=f"{os.path.dirname(xvg_file_paths[0])}/out_rms_mean.png")

if block and GUI:
    plt.show()