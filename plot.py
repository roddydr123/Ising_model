import matplotlib.pyplot as plt
import numpy as np
from ising_model import run_sim, updateSpinsGlauber, updateSpinsKawasaki


def single_magnetisation():
    data = np.genfromtxt('exam_data/ising.dat', delimiter=',', skip_header=1)
    plt.title("Magnetisation over time")
    plt.plot(data[:,0], data[:,1])
    plt.xlabel("Number of sweeps")
    plt.ylabel("Total magnetisation")
    plt.show()


def get_all_data():
    kT_list = np.arange(1, 3, 0.3)

    J = 1.0

    grid_size = 50
    # grid = np.ones([grid_size, grid_size])
    grid = np.random.choice([1, -1], [grid_size, grid_size])

    for kT in kT_list:
        filename = f"exam_data/temp-{kT}-K.dat"
        grid = run_sim(kT, updateSpinsKawasaki, 1000, grid_size, False, grid, J, filename)


def analyse():

    kT_list = np.arange(1, 3, 0.3)
    grid_size = 50

    fig = plt.figure(figsize=(12, 9))

    ax1 = fig.add_subplot(221)
    ax2 = fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4 = fig.add_subplot(224)

    susceptibilities = []
    mean_energies = []
    mean_magnetisations = []
    Cs = []

    for kT in kT_list:
        filename = f"exam_data/temp-{kT}-K.dat"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        sus = 1/(kT * grid_size**2) * np.var(data[:,1])
        susceptibilities.append(sus)

        C = 1/(kT**2 * grid_size**2) * np.var(data[:,2])
        Cs.append(C)

        mean_energies.append(np.average(data[:,2]))
        mean_magnetisations.append(np.abs(np.average(data[:,1])))

        # c_error = get_bootstrap_error(200, 500, data[:,2], grid_size**2 * kT)
        # x_error = get_bootstrap_error(200, 500, data[:,1], grid_size**2 * kT**2)
        c_error = get_jacknife_error(data[:,2], grid_size**2 * kT)
        x_error = get_jacknife_error(data[:,1], grid_size**2 * kT**2)
        e_error = np.std(np.abs(data[:,2]))
        m_error = np.std(np.abs(data[:,1]))



    ax1.errorbar(kT_list, mean_energies, yerr=e_error)
    ax1.set_ylabel("mean energy")
    ax2.errorbar(kT_list, mean_magnetisations, yerr=m_error)
    ax2.set_ylabel("mean magnetisation")
    ax3.errorbar(kT_list, Cs, yerr=c_error)
    ax3.set_ylabel("heat capacity")
    ax4.errorbar(kT_list, susceptibilities, yerr=x_error)
    ax4.set_ylabel("susceptibility")
    plt.show()



def get_bootstrap_error(n, k, data, constant):

    quantity = []

    for i in range(k):
        resample = np.random.choice(data, n)
        quantity.append(1/constant * np.var(resample))

    error = np.var(quantity)
    return error


def get_jacknife_error(data, constant) -> float:

    c = 1/constant * np.var(data)
    quantity = []

    for i in range(len(data)):
        resample = np.var(data[np.arange(len(data)) != i])
        quantity.append(1/constant * resample)

    error = np.sqrt(np.sum((quantity - c)**2))
    return error


get_all_data()
# analyse()