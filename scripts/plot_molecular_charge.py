import matplotlib.pyplot as plt
import numpy as np
from math import log
import os.path as path
import sys
# import glob
import csv
import time
import subprocess

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)



class Plotter:
    def __init__(self, mol):
        self.p = path.abspath(path.join(__file__ ,"../../"))
        # Inputs
        molfile = self.p+"/input/"+mol+".mol"
        self.mol = {'name': mol, 'file': molfile, 'mtime': path.getmtime(molfile)}
        self.atomlist = []

        # Outputs
        self.outDir = self.p + "/output/__Molecular/"+mol

        self.freeFile = self.outDir+"/freeDist.csv"
        self.boundFiles =[]
        self.intFile = self.outDir + "/intensity.csv"

        self.get_atoms()
        self.autorun=False

    def get_atoms(self):
        self.atomlist = []
        with open(self.mol['file'], 'r') as f:
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
                        self.atomlist.append({'name': a, 'file': file, 'mtime': path.getmtime(file)})
        self.boundFiles =[]
        for atomic in self.atomlist:
            self.boundFiles.append(self.outDir+"/dist_%s.csv"%atomic['name']);

    def update_inputs(self):
        self.get_atoms()
        self.mol['mtime'] = path.getmtime(self.mol['file'])

    def rerun_ac4dc(self):
        cmd = self.p+'/bin/rates '+self.mol['file']
        print("Running: ", cmd)
        subprocess.run(cmd, shell=True, check=True)

    def check_current(self):
        # Pull most recent atom mod time
        self.update_inputs()
        # Outputs are made at roughly the same time: Take the oldest.
        out_time = path.getmtime(self.freeFile)
        for bfile in self.boundFiles:
            out_time = min(out_time, path.getmtime(bfile))

        if self.mol['mtime'] > out_time:
            print("Master file %s is newer than most recent run" % self.mol['file'])
            return False

        for atomic in self.atomlist:
            if atomic['mtime'] > out_time:
                print("Dependent file %s is newer than most recent run" % atomic['file'])
                return False
        return True

    def update(self):
        raw = np.genfromtxt(self.intFile, comments='#')
        self.intensity = raw[:,1]
        self.time = raw[:, 0]
        with open(self.freeFile) as f:
            r = csv.reader(f, delimiter=' ', comment='#')
            erow = []
            for row in r:
                print(row)
                if row[0]=='#' and row[1]=='|':


        raw = np.genfromtxt(self.freeFile, comments='#')



    def setup_axes(self):
        fig.clf()
        ax1 = fig.add_subplot(111)
        ax2 = ax1.twinx()
        ax2.plot(self.time, self.intensity, lw = 2, c = 'black', ls = '--', alpha = 0.7)
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
