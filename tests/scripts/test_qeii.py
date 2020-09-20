import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import sys
import csv

raw = np.genfromtxt(sys.argv[1],comments='#',dtype=np.float64)
plt.plot(raw[0],raw[1])
plt.show()
