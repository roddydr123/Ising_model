import sys
import matplotlib.pyplot as plt
from physics import Ising

J=1.0
nstep=5000

if(len(sys.argv) != 4):
    print ("Usage python vis.py N T method")
    sys.exit()

lx=int(sys.argv[1]) 
ly=lx
kT=float(sys.argv[2])
method = sys.argv[3]

model = Ising(lx, kT, method, nstep, spins="random")
fig, ax = plt.subplots()
im=ax.imshow(model.spins, animated=True)

for n in range(nstep):
    for i in range(lx):
        for j in range(ly):
            if method == "G":
                model.updateSpinsGlauber()
            if method == "K":
                model.updateSpinsKawasaki()
            model.count += 1

    if(n%10==0):
        fig.canvas.flush_events()
        im=ax.imshow(model.spins, animated=True)
        ax.draw_artist(im)
        plt.pause(0.011)
