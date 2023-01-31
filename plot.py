import matplotlib.pyplot as plt
import numpy as np


read_data = np.genfromtxt("data/bs_errors.dat", delimiter=",")

plt.scatter(read_data[:,0], read_data[:,5])
# plt.errorbar(read_data[:,0], read_data[:,1], yerr=read_data[:,2])

plt.show()