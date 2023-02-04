import matplotlib.pyplot as plt
import numpy as np


read_data = np.genfromtxt("data/glauber.dat", delimiter=",")

fig = plt.figure(figsize=(12, 9))

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

ax1.scatter(read_data[:,0], (read_data[:,3]))
ax1.plot(read_data[:,0], (read_data[:,3]))
ax1.set_ylabel("Energy")
ax1.set_ylim(ax1.get_ylim()[0], 0)

ax2.scatter(read_data[:,0], (read_data[:,4]))
ax2.plot(read_data[:,0], (read_data[:,4]))
ax2.set_ylabel("Magnetisation")
ax2.set_ylim(0, ax2.get_ylim()[1])

ax3.scatter(read_data[:,0], (read_data[:,1]))
ax3.plot(read_data[:,0], (read_data[:,1]))
ax3.set_ylabel("Heat capacity")
ax3.set_ylim(0, ax3.get_ylim()[1])

ax4.scatter(read_data[:,0], (read_data[:,5]))
ax4.plot(read_data[:,0], (read_data[:,5]))
ax4.set_ylabel("Susceptibility")
ax4.set_ylim(0, ax4.get_ylim()[1])

plt.show()