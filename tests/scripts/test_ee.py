import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import sys


raw = np.genfromtxt(sys.argv[1])

#  colormap to use
cm = 'magma'

headers = []
with open(sys.argv[1], 'r') as f:
    while True:
        row = f.readline()
        print(row)
        if row[0] != '#':
            break
        else:
            headers.append(row)

Y = np.array(headers[-1].split(' ')[4:-1],dtype=np.float64)

t=raw[:,0]
dat = raw[:,1:]
# Y = np.arange(0, dat.shape[1])

def plot(min=0, max=None, log=False):
    Z = dat.T
    Z = np.ma.masked_where(Z < min, Z)
    if max is not None:
        Z = np.ma.masked_where(Z > max, Z)
    if log:
        Z = np.log(Z)

    fig, (ax) = plt.subplots(nrows=1,ncols=1)
    im = ax.pcolormesh(t, Y, Z, cmap = cm)
    ax.set_xlabel('Time')
    ax.set_ylabel('Energy')
    fig.colorbar(im)

def plot_sums(min=0, max=None):
    sumZ = np.sum(dat, axis=1)
    if max is None:
        max = np.max(sumZ)
    _fig, (ax) = plt.subplots(nrows=1,ncols=1)
    ax.plot(t, sumZ)
    ax.set_xlabel('Time')
    ax.set_ylabel('Total density')
    ax.set_ylim([0,max])


def maxwell(e, kT, n):
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def fit_step(idx = -1, guess = [1000, 1]):
    state = dat[idx,:]
    popt, _pcov = curve_fit(maxwell, Y, state, p0 = guess)
    print("kT = %f, n = %f" % tuple(popt))
    plt.plot(Y, state,label='t = %d' % t[idx])
    plt.plot(Y, maxwell(Y, popt[0], popt[1]), label='fit')

def plot_step(idx = -1, guess = [1000, 1]):
    state = dat[idx,:]
    plt.plot(Y, state,label='t = %d' % t[idx])

def plotsteps():
    for i in range(0,t.size,t.size//5):
        plot_step(i)
    plt.legend()