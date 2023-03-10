import matplotlib.pyplot as plt
import numpy as np


read_data = np.genfromtxt("findat/results.G3.dat", delimiter=",").T

fig = plt.figure(figsize=(12, 9))

ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)

# ax1.scatter(read_data[0], (read_data[3]))
ax1.errorbar(read_data[0], (read_data[3]), yerr=read_data[4], fmt="k-s", capsize=3)
ax1.set_ylabel("Energy")
ax1.set_xlabel("kT")
# ax1.set_ylim(ax1.get_ylim()[0], 0)

ax2.errorbar(read_data[0], (read_data[5]), yerr=read_data[6], fmt="k-s", capsize=3)
# ax2.plot(read_data[0], (read_data[4]))
ax2.set_ylabel("Magnetisation")
ax2.set_xlabel("kT")
# ax2.set_ylim(0, ax2.get_ylim()[1])

ax3.errorbar(read_data[0], (read_data[1]), yerr=read_data[2], fmt="k-s", capsize=3)
ax3.set_ylabel("Heat capacity")
ax3.set_xlabel("kT")
# ax3.set_ylim(0, ax3.get_ylim()[1])

ax4.errorbar(read_data[0], (read_data[7]), yerr=read_data[8], fmt="k-s", capsize=3)
# ax4.plot(read_data[0], (read_data[5]))
ax4.set_ylabel("Susceptibility")
ax4.set_xlabel("kT")
# ax4.set_ylim(0, ax4.get_ylim()[1])

plt.show()