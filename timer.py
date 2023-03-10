from physics import Ising
import numpy as np
from tqdm import tqdm
import timeit
import sys

J=1.0
nstep=10100
wait=100

lx=50
ly=lx

method = "G"

run_name = "g1"

def onerun(kT):

    model = Ising(lx, kT, method, nstep)

    # time = timeit.timeit(model.get_total_energy, number=100)
    en = model.get_total_energy()
    print(en)
    # print(time)
    sys.exit()

    for n in tqdm(range(nstep)):
        for i in range(lx):
            for j in range(ly):
                if method == "G":
                    model.updateSpinsGlauber()
                if method == "K":
                    model.updateSpinsKawasaki()
                model.count += 1

        if n % 10 == 0 and n > wait:

            model.get_total_magnetisation()
            model.get_total_energy()

    # with open(f"data/{model.temperature}.mags.{run_name}.dat", "w") as outfile:
    #     outfile.writelines(f"{model.magnetisation_list}")
    # with open(f"data/{model.temperature}.energies.{run_name}.dat", "w") as outfile:
    #     outfile.writelines(f"{model.energy_list}")

    sp_heat = model.get_heat_capacity()
    av_energy = np.average(model.energy_list)
    av_magnetisation = abs(np.average(model.magnetisation_list))
    susceptibility = model.get_susceptibility()

    c_error = model.get_bootstrap_error(200, 50, "c")
    x_error = model.get_bootstrap_error(200, 50, "x")

    # with open(f"data/results.{run_name}.dat", "a") as outfile:
    #     outfile.write(f"{kT}, {sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {susceptibility}, {x_error}\n")


def main():

    # clear file
    # f = open("data/results.dat", "w")
    # f.close()

    # make array of all the temperatures to loop through
    temps = np.round(np.arange(1, 3, 0.1), 3)
    # temps = [1.8]

    # for i in map(onerun, temps):
    #     pass

    for temp in temps:
        onerun(temp)


main()
