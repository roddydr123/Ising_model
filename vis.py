import sys
import matplotlib.pyplot as plt
from physics import Ising
import numpy as np

J=1.0
nstep=500

if(len(sys.argv) != 4):
    print ("Usage python ising.animation.py N T method")
    sys.exit()

lx=int(sys.argv[1]) 
ly=lx 
kT=float(sys.argv[2])
method = sys.argv[3]

model = Ising(lx, kT, method)
spins = model.spins
fig, ax = plt.subplots()
im=ax.imshow(spins, animated=True)

for n in range(nstep):
    for i in range(lx):
        for j in range(ly):
            if method == "G":
                model.updateSpinsGlauber()
            if method == "K":
                model.updateSpinsKawasaki()

    if(n%10==0): 
        fig.canvas.flush_events()
        spins = model.spins
        im=ax.imshow(spins, animated=True)
        ax.draw_artist(im)
        plt.pause(0.011)
        model.get_total_magnetisation()
        model.get_total_energy()

av_sp_heat = model.get_heat_capacity()
av_energy = np.average(model.energy_list)
av_magnetisation = np.average(model.magnetisation_list)
av_susceptibility = model.get_susceptibility()

c_error = model.get_bootstrap_error(200, 50, "c")

print(f"{kT}, {av_sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {av_susceptibility}")