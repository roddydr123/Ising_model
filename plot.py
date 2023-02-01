import matplotlib.pyplot as plt
import numpy as np


read_data = np.genfromtxt("data/results.dat", delimiter=",")

plt.scatter(read_data[:,0], np.abs(read_data[:,1]))
# plt.errorbar(read_data[:,0], read_data[:,1], yerr=read_data[:,2])

plt.show()