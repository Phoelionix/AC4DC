import numpy as np
import matplotlib.pyplot as plt
import sys


raw = np.genfromtxt(sys.argv[1])

#  colormap to use
cm = 'magma'

t=raw[:,0]
dat = raw[:,1:]
Y = np.arange(0, dat.shape[1])

fig, (ax) = plt.subplots(nrows=1,ncols=1)
im = ax.pcolormesh(t, Y, dat.T, cmap = cm)
ax.set_xlabel('Time')
ax.set_ylabel('Energy')
fig.colorbar(im)