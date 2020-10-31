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

from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)




# usage: python3 plot_slice_video.py [run handle] [num frames] [num initial]

plt.rcParams["font.size"] = 16

label = sys.argv[1] 

pl = Plotter(sys.argv[1])
pl.fig_steps.set_size_inches(8,5)

max_t = pl.timeData[-1]
min_t = pl.timeData[0]
nsteps = int(sys.argv[2])
times = np.linspace(min_t,max_t,nsteps)

pl.fig_steps.subplots_adjust(bottom=0.15,left=0.15,right=0.95,top=0.9)


vid_dir = 'video'

if not os.path.exists(vid_dir):
    os.makedirs(vid_dir)

n=0

num_init = 0
if len(sys.argv) > 3:
    num_init = int(sys.argv[3])
    n = num_init

# pad start with zeros to build suspense
before_times = times[:num_init + 1]
before_times = before_times-before_times[-1] + min_t
before_times = before_times[:-1]


# Formatting 
pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False
pl.ax_steps.set_ylim([1e-3, 1])
pl.ax_steps.set_xlim([1,10000]) 
pl.ax_steps.set_title('Carbon, 2mJ pulse, 100nm focus')
pl.ax_steps.loglog()

text = pl.ax_steps.text(2, 8e-2,'t = {:.2f} fs'.format(-10))

# Create a set of inset Axes: these should fill the bounding box allocated to
# them.
ax2 = pl.fig_steps.add_axes([0,0,1,1])
# Manually set the position and relative size of the inset axes within ax1
ip = InsetPosition(pl.ax_steps, [0.05,0.75,0.3,0.2])
ax2.set_axes_locator(ip)

pulseX = np.concatenate((before_times,pl.timeData))
pulseY = np.concatenate((np.zeros_like(before_times), pl.intensityData))

ax2.plot(pulseX, pulseY, 'k-')
ax2.yaxis.set_visible(False)
ax2.xaxis.set_visible(False)

pl.plot_fit(max_t, 2000, normed=False, lw=0.5, color='k')
pl.ax_steps.legend(loc='center left')
pl.plot_step(min_t, normed=False, lw=1, color='b')

for i in range(num_init):
    t = before_times[i]
    print(t)
    lineal = ax2.axvline(x=t,color='r')
    text.set_text('t = {:.2f} fs'.format(t))
    text.set_visible(True)
    plt.savefig('video/{:0>5}.png'.format(i))
    lineal.remove()

for t in times:
    print(t)
    vert = ax2.axvline(x=t,color='r')
    lines = pl.plot_step(t, normed=False, lw=1, color='b')
    text.set_text('t = {:.2f} fs'.format(t))
    plt.savefig('video/{:0>5}.png'.format(n))

    lines[-1].remove()
    vert.remove()
    
    n += 1

# pl.ax_steps.cla()
# pl.ax_steps.xaxis.get_major_formatter().labelOnlyBase = False
# pl.ax_steps.yaxis.get_major_formatter().labelOnlyBase = False
# pl.ax_steps.set_ylim([1e-3, 1])
# pl.ax_steps.set_xlim([1,10000])
    

    


