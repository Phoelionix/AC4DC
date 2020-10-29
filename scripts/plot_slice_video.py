import matplotlib
# matplotlib.use("pgf")
# matplotlib.rcParams.update({
#     "pgf.texsystem": "pdflatex",
#     'font.family': 'serif',
#     'text.usetex': True,
#     'pgf.rcfonts': False,
# })
import matplotlib.pyplot as plt
from plot_molecular_charge import Plotter
import sys
import numpy as np
import os

# usage: python3 plot_slice_video.py [run handle] [num frames]

plt.rcParams["font.size"] = 8

label = sys.argv[1] 

pl = Plotter(sys.argv[1])
pl.fig_steps.set_size_inches(5.3,3)

times = np.linspace(-10,20,int(sys.argv[2]))

pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.set_ylim([1e-3, 10])
pl.ax_steps.set_xlim([1,10000]) 
pl.fig_steps.subplots_adjust(bottom=0.15,left=0.2,right=0.95,top=0.2)


vid_dir = 'video'

if not os.path.exists(vid_dir):
    os.makedirs(vid_dir)

for n in range(len(times)):
    t = times[n]
    print(t)
    res = pl.plot_step(t, normed=False, lw=0.5)
    pl.ax_steps.set_title('t = {:.2} fs'.format(t))
    plt.savefig('video/{:0>5}.png'.format(n))
    res[-1].remove()




