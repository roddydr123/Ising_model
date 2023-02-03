from physics import Ising
import numpy as np
from tqdm import tqdm

J=1.0
nstep=10100
wait=100

lx=50
ly=lx

# n_randints = nstep * lx * ly

method = "G"

def onerun(kT):

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

    with open("energyfile.dat", "a") as file:
        file.write(model.energy_list)
    with open("magfile.dat", "a") as file:
        file.write(model.magnetisation_list)

    sp_heat = model.get_heat_capacity()
    av_energy = np.average(model.energy_list)
    av_magnetisation = abs(np.average(model.magnetisation_list))
    susceptibility = model.get_susceptibility()

    c_error = model.get_bootstrap_error(200, 50, "c")

    # print(f"{kT}, {sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {susceptibility}")
    with open("data/results.dat", "a") as outfile:
        outfile.write(f"{kT}, {sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {susceptibility}\n")


def main():

    # clear file
    f = open("data/results.dat", "w")
    f.close()

    # make array of all the temperatures to loop through
    temps = np.round(np.arange(1, 3, 0.3), 3)

    for temp in temps:
        onerun(temp)


main()