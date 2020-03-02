import matplotlib.pyplot as plt
import numpy as np
from math import log
import os.path as path
import sys
import time
import subprocess

plt.rcParams["font.size"] = 14
fig = plt.figure(figsize=(9, 6))

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

class Plotter:
    def __init__(self, species):
        self.p = path.abspath(path.join(__file__ ,"../../"))

        self.charge=np.zeros((1,1))
        self.intensity=np.zeros((1,1))
        self.charge_file= self.p+"/output/Charge_"+ species+".txt"
        self.intensity_file= self.p + "/output/Intensity_"+species+".txt"
        self.species=species
        self.update()
        self.methods = [f for f in dir(self) if callable(getattr(self, f))]
        self.methods = [f for f in self.methods if not str(f).startswith('__')]
        self.autorun = True

    def check_current(self):
        inp = self.p+"/input/"+self.species+".inp"
        oup = self.p+"/output/log_"+self.species+".txt"

        lastedit = path.getmtime(inp)
        lastrun = path.getmtime(oup)

        if lastrun < lastedit:
            print("File\n"+inp+"\n is newer than  \n"+oup)
            print("Rerunning ac4dc...")

            subprocess.run(self.p+'/bin/ac4dc input/{}.inp'.format(self.species), shell=True, check=True)

    def update(self):
        self.check_current()

        File = open(self.charge_file)
        # read the content into a list "Rates"
        charge = []
        for line in File:
            charge.append(line.split())
            for i, elem in enumerate(charge[-1]):
                charge[-1][i] = float(elem)
        File.close()
        charge = charge[100:-100]
        charge = np.array(charge)
        self.charge = charge.transpose()

        File = open(self.intensity_file)
        # read the content into a list "Rates"
        intensity = []
        for line in File:
            intensity.append(line.split())
            for i, elem in enumerate(intensity[-1]):
                intensity[-1][i] = float(elem)
        File.close()
        intensity = intensity[100:-100]
        intensity = np.array(intensity)
        self.intensity = intensity.transpose()

        self.Z = charge.shape[1] - 2
        # One for the timestamp, one for zero

    def plot_charges(self):
        fig.clf()
        ax = fig.add_subplot(111)

        for i, elem in enumerate(self.charge[:-1]):
            ax.plot(self.charge[-1], elem, lw = 2, alpha = 0.7, label = str(i))

        ax.plot(self.intensity[1], self.intensity[0], lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax.set_title("Charge state dynamics")
        ax.set_xlabel("Time, fs")
        ax.set_ylabel("Probability")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()

    def plot_avg_charge(self):
        fig.clf()
        ax = fig.add_subplot(111)
        expect = np.zeros(self.charge.shape[1])

        for i, elem in enumerate(self.charge[:-1]):
            expect += elem*i

        ax.plot(self.charge[-1], expect, lw = 2, alpha = 0.7, label = 'Average charge state')

        ax.plot(self.intensity[1], self.intensity[0], lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax.set_title("Charge state dynamics")
        ax.set_xlabel("Time, fs")
        ax.set_ylabel("Average charge")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()

    def plot_bound(self):
        fig.clf()
        ax = fig.add_subplot(111)
        expect = np.zeros(self.charge.shape[1])

        for i, elem in enumerate(self.charge[:-1]):
            expect += elem*i

        ax.plot(self.charge[-1], self.Z-expect, lw = 2, alpha = 0.7, label = 'Average charge state')

        ax.plot(self.intensity[1], self.intensity[0], lw = 2, c = 'black', ls = '--', alpha = 0.7)
        ax.set_title("Charge state dynamics")
        ax.set_xlabel("Time, fs")
        ax.set_ylabel("Average charge")

        plt.figlegend(loc = (0.11, 0.43))
        plt.subplots_adjust(left=0.1, right=0.92, top=0.93, bottom=0.1)

        plt.show()


pl = Plotter(sys.argv[1])
