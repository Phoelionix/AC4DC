import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import BSpline
from math import log
import os.path as path
import sys
# import glob
import csv
import time
import subprocess
import re

plt.style.use('seaborn')

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))

engine = re.compile(r'\d[spdf]\^\{(\d+)\}')

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def parse_elecs_from_latex(latexlike):
    # Parses a chemistry-style configuration to return total number of electrons
    # e.g. $1s^{1}2s^{1}$
    # ignores anything non-numerical or spdf
    num = 0
    for match in engine.finditer(latexlike):
        num += int(match.group(1))
    return num


ATOMS = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr'.split()
ATOMNO = {}
i = 1
for symbol in ATOMS:
    ATOMNO[symbol] = i
    i += 1
ATOMNO['forbidden_O'] = 8

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
                charge = ATOMNO[a] - parse_elecs_from_latex(states[i])
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
        fig.clf()
        ax = fig.add_subplot(111)
        ax2 = ax.twinx()
        ax2.plot(self.timeData, self.intensityData, lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax2.set_ylabel('Pulse fluence')
        ax.set_xlabel("Time, fs")
        return (ax, ax2)

    def plot_atom_total(self, a):
        ax, _ax2 = self.setup_axes()
        tot = np.sum(self.boundData[a], axis=1)
        ax.plot(self.timeData, tot)
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        plt.show()

    def plot_atom_raw(self, a):
        ax, _ax2 = self.setup_axes()
        for i in range(self.boundData[a].shape[1]):
            ax.plot(self.timeData, self.boundData[a][:,i], label = self.statedict[a][i])
        ax.set_title("Configurational dynamics")
        ax.set_ylabel("Density")
        fig.figlegend(loc = (0.11, 0.43))
        fig.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        fig.show()

    def plot_ffactor(self, a, num_tsteps = 20):
        fdists = np.genfromtxt('output/'+a+'/Xsections/Form_Factor.txt')
        # These correspond to the meaning of the FormFactor.txt entries themselves
        KMIN = 0
        KMAX = 2
        dim = len(fdists.shape)
        kgrid = np.linspace(KMIN,KMAX,fdists.shape[0 if dim == 1 else 1])
        fig2 = plt.figure()
        ax = fig2.add_subplot(111)
        ax.set_xlabel('k, atomic units')
        
        timedata = self.boundData[a][:,:-1] # -1 excludes the bare nucleus
        dynamic_k = np.tensordot(fdists.T, timedata.T,axes=1) 
        step = dynamic_k.shape[1] // num_tsteps
        for i in range(0, dynamic_k.shape[1], step):
            ax.plot(kgrid, dynamic_k[:,i])
        fig.show()

    def plot_charges(self, a):
        ax, _ax2 = self.setup_axes()
        self.aggregate_charges()
        for i in range(self.chargeData[a].shape[1]):
            max_at_zero = np.max(self.chargeData[a][0,:])
            mask = self.chargeData[a][:,i] > max_at_zero*2
            mask |= self.chargeData[a][:,i] < -max_at_zero*2
            Y = np.ma.masked_where(mask, self.chargeData[a][:,i])
            ax.plot(self.timeData, Y, label = "%d+" % i)
        ax.set_title("Charge state dynamics")
        ax.set_ylabel("Proportion (arb. units)")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        plt.show()

    def plot_tot_charge(self):
        ax, _ax2 = self.setup_axes()
        self.aggregate_charges()
        Q = np.zeros(self.timeData.shape[0])
        for a in self.atomdict:
            atomic_charge = np.zeros(self.timeData.shape[0])
            for i in range(self.chargeData[a].shape[1]):
                atomic_charge += self.chargeData[a][:,i]*i
            ax.plot(self.timeData, atomic_charge, label = a)
            Q += atomic_charge

        tot_free_Q =-1*np.sum(self.freeData, axis=1)*(self.energyKnot[6]-self.energyKnot[5])
        ax.plot(self.timeData, tot_free_Q, label = 'Free')
        print(tot_free_Q/Q)
        Q += tot_free_Q
        ax.set_title("Charge state dynamics")
        ax.plot(self.timeData, Q, label='total')
        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)
        plt.show()


    def plot_all_charges(self):
        for a in self.atomdict:
            self.plot_charges(a)

    def plot_free(self, N=100, log=False, min = None, max=None, tmax = None):
        fig.clf()
        # Need to turn freeData (matrix of BSpline coeffs)
        # freeeData [t, c]
        # into matrix of values.
        if tmax is None:
            Z = self.freeData.T
            T = self.timeData
        else:
            Z = self.freeData.T [:, :tmax]
            T = self.timeData [:tmax]
        
        if log:
            Z = np.log(Z)
        
        if min is not None:
            Z = np.ma.masked_where(Z < min, Z)

        if max is not None:
            Z = np.ma.masked_where(Z > max, Z)

        plt.contourf(T, self.energyKnot, Z, N, cmap='viridis')
        plt.title("Free electron energy distribution")
        plt.ylabel("Energy (eV)")
        plt.xlabel("Time, fs")
        plt.show()
        plt.colorbar()




pl = Plotter(sys.argv[1])
pl.plot_tot_charge()
