import matplotlib.rcsetup as rcsetup
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BSpline
from math import log
import os.path as path
import matplotlib.colors as colors
import sys
# import glob
import csv
import time
import subprocess
import re
from matplotlib.ticker import LogFormatter 
import random
from scipy.optimize import curve_fit
from scipy.stats import linregress


plt.style.use('seaborn-muted')
plt.rcParams["font.size"] = 9

engine = re.compile(r'(\d[spdf])\^\{(\d+)\}')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_elecs_from_latex(latexlike):
    # Parses a chemistry-style configuration to return total number of electrons
    # e.g. $1s^{1}2s^{1}$
    # ignores anything non-numerical or spdf
    qdict = {}
    for match in engine.finditer(latexlike):
        qdict[match.group(1)] = int(match.group(2))
    return qdict


ATOMS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split()
ATOMNO = {}
i = 1
for symbol in ATOMS:
    ATOMNO[symbol] = i
    i += 1
ATOMNO['forbidden_O'] = 8

def get_colors(num, seed):
    idx = list(np.linspace(0, 1, num+1))[1:]
    random.seed(seed)
    random.shuffle(idx)
    idx.insert(0,0)
    C = plt.cm.nipy_spectral
    return C(idx)

class Plotter:
    # Initialistation: Plotter(water)
    # expects there to be a control file named input/water.mol
    def __init__(self, mol):
        self.p = path.abspath(path.join(__file__ ,"../../"))
        # Inputs
        molfile = self.p+"/input/"+mol+".mol"
        self.mol = {'name': mol, 'infile': molfile, 'mtime': path.getmtime(molfile)}
        # Stores the atomic input files read by ac4dc
        self.atomdict = {}
        self.statedict = {}

        # Outputs
        self.outDir = self.p + "/output/__Molecular/" + mol
        self.freeFile = self.outDir+"/freeDist.csv"
        self.intFile = self.outDir + "/intensity.csv"

        self.boundData={}
        self.chargeData={}
        self.freeData=None
        self.intensityData=None
        self.energyKnot=None
        self.timeData=None

        self.get_atoms()
        self.update_outputs()
        self.autorun=False

        self.setup_step_axes()


    def setup_step_axes(self):
        self.fig_steps = plt.figure(figsize=(4,3))
        self.ax_steps = self.fig_steps.add_subplot(111)
        self.fig_steps.subplots_adjust(left=0.18,right=0.97, top=0.97, bottom= 0.16)
        

    # Reads the control file specified by self.mol['infile']
    # and populates the atomdict data structure accordingly
    def get_atoms(self):
        self.atomdict = {}
        with open(self.mol['infile'], 'r') as f:
            reading = False
            for line in f:
                if line.startswith("#ATOMS"):
                    reading=True
                    continue
                elif line.startswith("#"):
                    reading=False
                    continue
                if reading:
                    a = line.split(' ')[0].strip()
                    if len(a) != 0:
                        file = self.p + '/input/' + a + '.inp'
                        self.atomdict[a]={
                            'infile': file,
                            'mtime': path.getmtime(file),
                            'outfile': self.outDir+"/dist_%s.csv"%a}

    def update_inputs(self):
        self.get_atoms()
        self.mol['mtime'] = path.getmtime(self.mol['infile'])

    def rerun_ac4dc(self):
        cmd = self.p+'/bin/ac4dc2 '+self.mol['infile']
        print("Running: ", cmd)
        subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        # Pull most recent atom mod time
        self.update_inputs()
        # Outputs are made at roughly the same time: Take the oldest.
        out_time = path.getmtime(self.freeFile)
        for atomic in self.atomdict.values():
            out_time = min(out_time, path.getmtime(atomic['outfile']))

        if self.mol['mtime'] > out_time:
            print("Master file %s is newer than most recent run" % self.mol['infile'])
            return False

        for atomic in self.atomdict.values():
            if atomic['mtime'] > out_time:
                print("Dependent input file %s is newer than most recent run" % atomic['infile'])
                return False
        return True

    def get_free_energy_spec(self):
        erow = []
        with open(self.freeFile) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find energy grid specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        erow.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return erow

    def get_bound_config_spec(self, a):
        # gets the configuration strings corresponding to atom a
        specs = []
        with open(self.atomdict[a]['outfile']) as f:
            r = csv.reader(f, delimiter=' ')
            for row in r:
                if row[0] != '#':
                    raise Exception('parser could not find atomic state specification, expected #  |')
                reading=False
                for entry in row:
                    # skip whitespace
                    if entry == '' or entry == '#':
                        continue
                    # it's gotta be a pipe or nothin'
                    if reading:
                        specs.append(entry)
                    elif entry == '|':
                        reading=True
                    else:
                        break
                if reading:
                    break
        return specs

    def aggregate_charges(self):
        # populates self.chargeData based on contents of self.boundData
        for a in self.atomdict:
            states = self.statedict[a]
            if len(states) != self.boundData[a].shape[1]:
                msg = 'states parsed from file header disagrees with width of data width'
                msg += ' (got %d, expected %d)' % (len(states), self.boundData[a].shape[1])
                raise RuntimeError(msg)

            self.chargeData[a] = np.zeros((self.boundData[a].shape[0], ATOMNO[a]+1))
            for i in range(len(states)):
                orboccs = parse_elecs_from_latex(states[i])
                charge = ATOMNO[a] - sum(orboccs.values())
                self.chargeData[a][:, charge] += self.boundData[a][:, i]


    def update_outputs(self):
        raw = np.genfromtxt(self.intFile, comments='#', dtype=np.float64)
        self.intensityData = raw[:,1]
        self.timeData = raw[:, 0]
        self.energyKnot = np.array(self.get_free_energy_spec(), dtype=np.float64)
        raw = np.genfromtxt(self.freeFile, comments='#', dtype=np.float64)
        self.freeData = raw[:,1:]
        for a in self.atomdict:
            raw = np.genfromtxt(self.atomdict[a]['outfile'], comments='#', dtype=np.float64)
            self.boundData[a] = raw[:, 1:]
            self.statedict[a] = self.get_bound_config_spec(a)

    def go(self):
        if not self.check_current():
            self.rerun_ac4dc()
            self.update_outputs()

    # makes a blank plot showing the intensity curve
    def setup_axes(self):
        self.fig = plt.figure(figsize=(3,2.5))
        ax = self.fig.add_subplot(111)
        ax2 = ax.twinx()
        ax2.plot(self.timeData, self.intensityData, lw = 1, c = 'black', ls = ':', alpha = 0.7)
        # ax2.set_ylabel('Pulse Intensity (photons cm$^{-2}$ s$^{-1}$)')
        ax.set_xlabel("Time (fs)")
        ax.tick_params(direction='in')
        ax2.axes.get_yaxis().set_visible(False)
        # ax2.tick_params(None)
        ax.get_xaxis().get_major_formatter().labelOnlyBase = False
        return (ax, ax2)

    def plot_atom_total(self, a):
        ax, _ax2 = self.setup_axes()
        tot = np.sum(self.boundData[a], axis=1)
        ax.plot(self.timeData, tot)
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        self.fig.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        
    def plot_atom_raw(self, a):
        ax, _ax2 = self.setup_axes()
        for i in range(self.boundData[a].shape[1]):
            ax.plot(self.timeData, self.boundData[a][:,i], label = self.statedict[a][i])
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        self.fig.subplots_adjust(left=0.2, right=0.92, top=0.93, bottom=0.1)

    def plot_ffactor(self, a, num_tsteps = 10):
        fdists = np.genfromtxt('output/'+a+'/Xsections/Form_Factor.txt')
        # These correspond to the meaning of the FormFactor.txt entries themselves
        KMIN = 0
        KMAX = 2
        dim = len(fdists.shape)
        kgrid = np.linspace(KMIN,KMAX,fdists.shape[0 if dim == 1 else 1])
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.set_xlabel('k, atomic units')
        ax.set_ylabel('Elastic scattering form factor (arbitrary units)')
        
        timedata = self.boundData[a][:,:-1] # -1 excludes the bare nucleus
        dynamic_k = np.tensordot(fdists.T, timedata.T,axes=1) 
        step = dynamic_k.shape[1] // num_tsteps
        for i in range(0, dynamic_k.shape[1], step):
            ax.plot(kgrid, dynamic_k[:,i])
        fig2.legend()
        fig2.show()

    def plot_charges(self, a, rseed=404):
        ax, _ax2 = self.setup_axes()
        self.aggregate_charges()
        
        ax.set_prop_cycle(rcsetup.cycler('color', get_colors(self.chargeData[a].shape[1],rseed)))
        for i in range(self.chargeData[a].shape[1]):
            max_at_zero = np.max(self.chargeData[a][0,:])
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])
            ax.plot(self.timeData, Y, label = "%d+" % i)
        # ax.set_title("Charge state dynamics")
        ax.set_ylabel(r"Density (\AA$^{-3}$)")

        self.fig.legend(loc = "right")
        self.fig.subplots_adjust(left=0.11, right=0.81, top=0.93, bottom=0.1)

    def plot_subshell(self, a, subshell='1s',rseed=404):
        if not hasattr(self, 'ax_subshell'):
            self.ax_subshell, _ax2 = self.setup_axes()

        states = self.statedict[a]
        tot = np.zeros_like(self.boundData[a][:,0])
        for i in range(len(states)):
            orboccs = parse_elecs_from_latex(states[i])
            tot += orboccs[subshell]*self.boundData[a][:,i]

        self.ax_subshell.plot(self.timeData, tot, label=subshell)
        # ax.set_title("Charge state dynamics")
        self.ax_subshell.set_ylabel(r"Electron density (\AA$^{-3}$)")

        self.fig.subplots_adjust(left=0.2, right=0.95, top=0.95, bottom=0.2)


    def plot_tot_charge(self, every=1):
        ax, _ax2 = self.setup_axes()
        self.aggregate_charges()
        self.fig.subplots_adjust(left=0.22, right=0.95, top=0.95, bottom=0.17)

        T = self.timeData[::every]
        self.Q = np.zeros(T.shape[0])
        for a in self.atomdict:
            atomic_charge = np.zeros(T.shape[0])
            for i in range(self.chargeData[a].shape[1]):
                atomic_charge += self.chargeData[a][::every,i]*i
            ax.plot(T, atomic_charge, label = a)
            self.Q += atomic_charge

        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        tot_free_Q =-1*np.dot(self.freeData, de)
        ax.plot(T, tot_free_Q[::every], label = 'Free')
        ax.set_ylabel("Charge density ($e$ \AA$^{-3}$)")
        self.Q += tot_free_Q[::every]
        # ax.set_title("Charge Conservation")
        ax.plot(T, self.Q, label='total')
        ax.legend(loc = 'upper left')
        return ax


    def plot_all_charges(self, rseed=404):
        for a in self.atomdict:
            self.plot_charges(a,rseed)

    def plot_free(self, N=100, log=False, min = 0, max=None, every = None):
        self.fig_free = plt.figure(figsize=(3.1,2.5))
        ax = self.fig_free.add_subplot(111)
        self.fig_free.subplots_adjust(left=0.12, top=0.96, bottom=0.16,right=0.95)
        # Need to turn freeData (matrix of BSpline coeffs)
        # freeeData [t, c]
        # into matrix of values.
        if every is None:
            Z = self.freeData.T
            T = self.timeData
        else:
            Z = self.freeData.T [:, ::every]
            T = self.timeData [::every]
        
        Z = np.ma.masked_where(Z < min, Z)

        if max is not None:
            Z = np.ma.masked_where(Z > max, Z)
        else:
            max = Z.max()

        # if log:
        #     Z = np.log10(Z)

        
        ax.set_facecolor('black')

        if log:
            cm = ax.pcolormesh(T, self.energyKnot*1e-3, Z, shading='gouraud',norm=colors.LogNorm(vmin=min, vmax=max),cmap='magma',rasterized=True)
            cbar = self.fig_free.colorbar(cm)
        else:
            cm = ax.contourf(T, self.energyKnot, Z, N, cmap='magma',rasterized=True)
            cbar = self.fig_free.colorbar(cm)

        ax.set_ylabel("Energy (keV)")
        ax.set_xlabel("Time (fs)")
        
        # if log:
        #     minval = np.floor(np.min(Z, axis=(0,1)))
        #     maxval = np.ceil(np.max(Z,axis=(0,1)))
        #     formatter = LogFormatter(10, labelOnlyBase=True) 
            
        #     cbar.ax.yaxis.set_ticks(vals)
        #     cbar.ax.yaxis.set_ticklabels(10**vals)
        # else:
        #     cbar = self.fig_free.colorbar(cm)
        cbar.ax.set_ylabel('Free Electron Density, Å$^{-3}$', rotation=270,labelpad=20)

        # plot the intensity
        ax2 = ax.twinx()
        ax2.plot(T, self.intensityData[::every],'w:',linewidth=1)

        ax2.get_yaxis().set_visible(False)

        

    def plot_free_raw(self, N=100, log=False, min = None, max=None):
        plt.figure()
        rawdata = np.genfromtxt(self.outDir+'/freeDistRaw.csv')
        T = rawdata[:,0]
        self.rawZ = rawdata[:,1:].T
        Z = self.rawZ

        if log:
            Z = np.log(Z)
        
        if min is not None:
            Z = np.ma.masked_where(Z < min, Z)

        if max is not None:
            Z = np.ma.masked_where(Z > max, Z)

        Yax = np.arange(Z.shape[0])
        plt.contourf(T, Yax, Z, N, shading='nearest',cmap='magma')
        plt.title("Free electron energy distribution")
        plt.ylabel("Energy (eV)")
        plt.xlabel("Time, fs")
        plt.show()
        plt.colorbar()

    def plot_step(self, t, normed=True, fitE=None, c=None):        
        self.ax_steps.set_xlabel('Energy, eV')
        self.ax_steps.set_ylabel('$f(\\epsilon) \\Delta \\epsilon$')
        self.ax_steps.loglog()
        # norm = np.sum(self.freeData[n,:])
        n = self.timeData.searchsorted(t)
        data = self.freeData[n,:]
        X = self.energyKnot

        if normed:
            tot = np.sum(data)
            data /= tot
        
        return self.ax_steps.plot(X, data*X, label='%1.1f fs' % t)

    def plot_fit(self, t, fitE, normed=True, c=None):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        data = self.freeData[t_idx,:]
        if normed:
            tot = np.sum(data)
            data /= tot

        Xdata = self.energyKnot[:fit]
        Ydata = data[:fit]
        mask = np.where(Ydata > 0)
        T, n = fit_maxwell(Xdata, Ydata)
        return self.ax_steps.plot(self.energyKnot, 
            maxwell(self.energyKnot, T, n)*self.energyKnot,
            '--', lw=1,label='%3.1f eV' % T, color=c)

    def get_temp(self, t, fitE):
        t_idx = self.timeData.searchsorted(t)
        fit = self.energyKnot.searchsorted(fitE)
        Xdata = self.energyKnot[:fit]
        Ydata = self.freeData[t_idx,:fit]
        T, n = fit_maxwell(Xdata, Ydata)
        return (T, n)

    def get_density(self, t):
        t_idx = self.timeData.searchsorted(t)
        de = np.append(self.energyKnot, self.energyKnot[-1]*2 - self.energyKnot[-2])
        de = de [1:] - de[:-1]
        return np.dot(self.freeData[t_idx, :], de)

        

def fit_maxwell(X, Y):
    guess = [200, 12]
    popt, _pcov = curve_fit(maxwell, X, Y, p0 = guess)
    return popt

def maxwell(e, kT, n):
    return n * np.sqrt(e/(np.pi*kT**3)) * np.exp(-e/kT)

def lnmaxwell(e, kT, n):
    return np.log(n) + 0.5*np.log(e/np.pi*kT**3) - e /kT

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

if __name__ == "__main__":
    pl = Plotter(sys.argv[1])
    pl.plot_free(log=True,min=1e-7)
    # pl.plot_all_charges()
    plt.show()
