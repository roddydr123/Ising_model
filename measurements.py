from physics import Ising
import numpy as np
from tqdm import tqdm

J=1.0
nstep=150
wait=100

lx=50
ly=lx

method = "G"

run_name = "test"

write_to_file = False

def onerun(kT, model):

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

    sp_heat = model.get_heat_capacity()
    av_energy = np.average(model.energy_list)
    av_magnetisation = abs(np.average(model.magnetisation_list))
    susceptibility = model.get_susceptibility()

    c_error = model.get_bootstrap_error(200, 50, "c")
    x_error = model.get_bootstrap_error(200, 50, "x")
    e_error = np.std(model.energy_list)
    m_error = np.std(model.magnetisation_list)

    if write_to_file is True:
        with open(f"data/{model.temperature}.mags.{run_name}.dat", "w") as outfile:
            outfile.writelines(f"{model.magnetisation_list}")
        with open(f"data/{model.temperature}.energies.{run_name}.dat", "w") as outfile:
            outfile.writelines(f"{model.energy_list}")
        with open(f"data/results.{run_name}.dat", "a") as outfile:
            outfile.write(f"{kT}, {sp_heat}, {c_error}, {av_energy}, {e_error}, {av_magnetisation}, {m_error}, {susceptibility}, {x_error}\n")
    return model.spins


def main():

    # make array of all the temperatures to loop through
    temps = np.round(np.arange(1, 3, 0.1), 3)

    for i, temp in enumerate(temps):
        # if its the first temperature, generate a model with new spins.
        # if its not the first, pass the spins matrix from the previous run.
        if i == 0:
            model = Ising(lx, temp, method, nstep)
        else:
            model = Ising(lx, temp, method, nstep, spins)
        spins = onerun(temp, model)
        


main()