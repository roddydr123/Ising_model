import sys
import matplotlib.pyplot as plt
from physics import Ising
import numpy as np

J=1.0
nstep=10000
filename = "output.dat"

if(len(sys.argv) != 4):
    print ("Usage python ising.animation.py N T method")
    sys.exit()

lx=int(sys.argv[1]) 
ly=lx 
kT=float(sys.argv[2])
method = sys.argv[3]

# clear outfile
with open(filename, "w") as f:
    f.write("magnetisation, energy\n")

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
                
#occasionally plot or update measurements, eg every 10 sweeps
    if(n%10==0): 
        # update measurements
        # dump output
        # f=open('spins.dat','w')
        # for i in range(lx):
        #     for j in range(ly):
        #         f.write('%d %d %lf\n'%(i,j,spins[i,j]))
        # f.close()
        fig.canvas.flush_events()
        spins = model.spins
        im=ax.imshow(spins, animated=True)
        ax.draw_artist(im)
        plt.pause(0.011)

        # calculate data for storage
        total_magnetisation = model.get_total_magnetisation()

        total_energy = model.get_total_energy()
        
        # record required data
        with open(filename, "a") as f:
            f.write(f"{total_magnetisation}, {total_energy}\n")
