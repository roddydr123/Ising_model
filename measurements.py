from physics import Ising
import numpy as np
from tqdm import tqdm

J=1.0
nstep=10100
wait=100

lx=50
ly=lx

method = "G"

def onerun(kT):
    # no visualisation
    # run a sim for 10000 sweeps
    # start measurements after 100 sweeps
    # record the magnetisation and susceptibility every 10 sweeps

    model = Ising(lx, kT, method)

    for n in tqdm(range(nstep)):
        for i in range(lx):
            for j in range(ly):
                if method == "G":
                    model.updateSpinsGlauber()
                if method == "K":
                    model.updateSpinsKawasaki()

        if n % 10 == 0 and n > wait:

            model.get_total_magnetisation()
            model.get_total_energy()

    av_sp_heat = model.get_heat_capacity()
    av_energy = np.average(model.energy_list)
    av_magnetisation = np.average(model.magnetisation_list)
    av_susceptibility = model.get_susceptibility()

    c_error = model.get_bootstrap_error(200, 50, "c")
    with open("data/results.dat", "a") as outfile:
        outfile.write(f"{kT}, {av_sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {av_susceptibility}\n")


def main():

    # clear file
    f = open("data/results.dat", "w")
    f.close()

    # make array of all the temperatures to loop through
    temps = np.round(np.arange(1, 3, 0.1), 3)

    for temp in tqdm(temps):
        onerun(temp)


main()