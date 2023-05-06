import matplotlib.pyplot as plt
import numpy as np


def magnetisation():
    data = np.genfromtxt('ising.dat', delimiter=',', skip_header=1)
    plt.title("Magnetisation over time")
    plt.plot(data[:,0], data[:,1])
    plt.xlabel("Number of sweeps")
    plt.ylabel("Total magnetisation")
    plt.show()

magnetisation()