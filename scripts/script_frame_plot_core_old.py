import stringprep
import matplotlib.rcsetup as rcsetup
import matplotlib.pyplot as plt
import numpy as np
# from scipy.interpolate import BSpline
from math import log
import os.path as path
import os
import matplotlib.colors as colors
import sys
# import glob
import csv
import subprocess
import re
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress

plt.rcParams["font.size"] = 9

def get_colors(num, seed):
    idx = list(np.linspace(0, 1, num))[1:]
    random.seed(seed)
    # random.shuffle(idx)
    idx.insert(0,0)
    C = plt.get_cmap('nipy_spectral')
    return C(idx)

class Plotter:
    # Example initialisation: Plotter(water)
    # --> expects there to be a control file named water.mol within AC4DC/input/ or a subdirectory.
    def __init__(self):
        self.setup_step_axes()
    def update(self,energies,densities):
        self.density = np.array(densities)
        self.energyKnot = np.array(energies)         

    def setup_step_axes(self):
        self.fig_steps = plt.figure(figsize=(4,3))
        self.ax_steps = self.fig_steps.add_subplot(111)
        self.fig_steps.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)

    # Plots a single point in time.
    def plot_step(self, normed=True, **kwargs):        
        self.ax_steps.set_xlabel('Energy (eV)')
        self.ax_steps.set_ylabel('$f(\\epsilon) \\Delta \\epsilon$')
        self.ax_steps.loglog()
        # norm = np.sum(self.freeData[n,:])
        dens = self.density
        E = self.energyKnot

        if normed:
            tot = self.get_density()
            dens /= tot
            dens /= 4*3.14
        
        return self.ax_steps.plot(E, dens*E, **kwargs) # label='%1.1f fs' % t

    def plot_fit(self, fitE, normed=True, **kwargs):
        fit = self.energyKnot.searchsorted(fitE)
        dens = self.density
        if normed:
            tot = self.get_density()
            dens /= tot
            dens/=4*3.14

        Xdata = self.energyKnot[:fit]
        Ydata = dens[:fit]
        T, n = fit_maxwell(Xdata, Ydata)
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, T, n)*self.energyKnot,
            '--',label='%3.1f eV' % T, **kwargs)

    def plot_maxwell(self, kT, n, **kwargs):
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, kT, n)*self.energyKnot,
            '--',label='%3.1f eV' % kT, **kwargs)

    def get_density(self):
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        return np.dot(self.density, de)

        

def fit_maxwell(X, Y):
    guess = [200, 12]
    # popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess, sigma=1/(X+10))
    popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
    return popt

def maxwell(e, kT, n):
    if kT < 0:
        return 0 # Dirty silencing of fitting error - note we get negative values from unphysical oscillations, so this increases the average value around this point. -S.P.
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def lnmaxwell(e, kT, n):
    return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
