import matplotlib.pyplot as plt
import numpy as np
from math import log
import os.path as path
import sys
import glob
import time
import subprocess

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_atom_files(filename):
    atomlist = []
    with open(filename, 'r') as f:
        reading = False
        for line in f:
            if line.startswith("#ATOMS"):
                reading=True
                continue
            elif line.startswith("#"):
                reading=False
                continue
            if reading:
                atomlist.append(line.split(' ')[0])

    retval = []
    for atom in atomlist:
        retval.append(("/input/"+atom+".inp", None))
    for pair in retval:
        pair[1] = path.getmtime(pair[0])
    return retval

class Plotter:
    def __init__(self, species):
        self.p = path.abspath(path.join(__file__ ,"../../"))
        dataDir = self.p + "output/__Molecular/"+species
        self.freeFile = dataDir+"/freeDist.csv"
        self.boundFiles =[]
        for file in glob.glob(self.p + "output/__Molecular/"+species+"/dist_*.csv"):
            self.boundFiles.append(file);
        self.autorun=False

    def rerun_ac4dc(self):
        cmd = self.p+'/bin/rates'+self.p+'input/{}.inp'.format(self.species)
        print("Running: ", cmd)
        subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        inp = self.p+"/input/"+self.species+".mol"

        lastedit = path.getmtime(inp)
        lastrun = path.getmtime(oup)

        if lastrun < lastedit:
            print("Master file '"+inp+"' is newer than '"+oup+"'")
            self.rerun_ac4dc()
            return

        # check if the atom files have been newly edited
        for name, ledit in get_atom_files(inp):
            if lastrun < ledit:
                print("Dependent file '"+name+"' is newer than '"+oup+"'")
                self.rerun_ac4dc()
                return

    def update(self):
        self.check_current()
        with open(self.charge_file) as f:


    def setup_axes(self):
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        ax2.plot(self.intensity[1], self.intensity[0], lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax2.set_ylabel('Pulse fluence')
        ax1.set_xlabel("Time, fs")
        return (ax1, ax2)

    def plot_charges(self):
        ax, pulseax = self.setup_axes()

        for i, elem in enumerate(self.charge[:-1]):
            ax.plot(self.charge[-1], elem, lw = 2, alpha = 0.7, label = str(i))

        ax.set_title("Charge state dynamics")
        ax.set_ylabel("Probability")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()


    def aggregate_charges(self):
        expect = np.zeros(self.charge.shape[1])

        for i, elem in enumerate(self.charge[:-1]):
            expect += elem*i
        return expect

    def plot_avg_charge(self):
        ax, pulseax = self.setup_axes()

        ax.plot(self.charge[-1], self.aggregate_charges(), lw = 2, alpha = 0.7, label = 'Average charge state')
        ax.set_title("Charge state dynamics")
        ax.set_ylabel("Average charge")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()

    def plot_bound(self):
        ax, pulseax = self.setup_axes()

        ax.plot(self.charge[-1], self.Z-self.aggregate_charges(), lw = 2, alpha = 0.7, label = 'Average charge state')
        ax.set_title("Charge state dynamics")
        ax.set_ylabel("Average charge")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()


pl = Plotter(sys.argv[1])
