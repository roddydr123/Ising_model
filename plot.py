import matplotlib.pyplot as plt
import numpy as np


def magnetisation():
    data = np.genfromtxt('ising.dat', delimiter=',', skip_header=1)
    plt.title("Magnetisation over time")
    plt.plot(data[:,0], data[:,1])
    plt.xlabel("Number of sweeps")
    plt.ylabel("Total magnetisation")
    plt.show()


    # sp_heat = get_heat_capacity()
    # av_energy = np.average(energy_list)
    # av_magnetisation = abs(np.average(magnetisation_list))
    # susceptibility = get_susceptibility()

    # c_error = get_bootstrap_error(200, 50, "c")
    # x_error = get_bootstrap_error(200, 50, "x")
    # e_error = np.std(np.abs(energy_list))
    # m_error = np.std(np.abs(magnetisation_list))

    # if write_to_file is True:
    #     with open(f"data/{model.temperature}.mags.{run_name}.dat", "w") as outfile:
    #         for line in model.magnetisation_list:
    #             outfile.write(str(line)+"\n")
    #     with open(f"data/{model.temperature}.energies.{run_name}.dat", "w") as outfile:
    #         for line in model.energy_list:
    #             outfile.write(str(line)+"\n")
    #     with open(f"data/results.{run_name}.dat", "a") as outfile:
    #         outfile.write(f"{kT}, {sp_heat}, {c_error}, {av_energy}, {e_error}, {av_magnetisation}, {m_error}, {susceptibility}, {x_error}\n")
    # return model.spins



magnetisation()