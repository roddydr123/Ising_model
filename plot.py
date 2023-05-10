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
    """Function to run simulation with multiple values of kT and save results."""

    kT_list = np.arange(1, 3, 0.3)

    J = 1.0

    grid_size = 50
    # grid = np.ones([grid_size, grid_size])
    grid = np.random.choice([1, -1], [grid_size, grid_size])

    for kT in kT_list:
        filename = f"exam_data/temp-{kT}-K.dat"
        grid = run_sim(kT, updateSpinsKawasaki, 10100, grid_size, False, grid, J, filename)


def analyse():
    """Load and plot the saved data from multiple kT simulations."""
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
        filename = f"exam_data/temp-{kT}.dat"
        data = np.genfromtxt(filename, delimiter=',', skip_header=1)
        sus = 1/(kT * grid_size**2) * np.var(data[:,1])
        susceptibilities.append(sus)

        C = 1/(kT**2 * grid_size**2) * np.var(data[:,2])
        Cs.append(C)

        mean_energies.append(np.average(data[:,2]))
        mean_magnetisations.append(np.abs(np.average(data[:,1])))

        c_error = get_jacknife_error(data[:,2], lambda x: np.var(x) / (grid_size**2 * kT**2))
        x_error = get_jacknife_error(data[:,1], lambda x: np.var(x) / (grid_size**2 * kT))
        # c_error = get_bootstrap_error(data[:,2], lambda x: np.var(x) / (grid_size**2 * kT**2))
        # x_error = get_bootstrap_error(data[:,1], lambda x: np.var(x) / (grid_size**2 * kT))
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


def get_jacknife_error(data, interesting_quantity_func):
    """
    Working jacknife error function.
    
    data: raw data from which interesting quantity is calculated.
    
    interesting quantity func: function which takes raw data and returns interesting quantity.
    """

    # make sure data to be resampled is an array.
    data = np.array(data)

    # calculate the interested quantity using all data.
    all_samples = interesting_quantity_func(data)

    # prepare an array for the quantities calculated with resampling.
    resampled_quantity = np.zeros(len(data))

    # loop over all the data and remove one sample each time.
    for i in range(len(data)):

        # array with all but one data point in it.
        resample = data[np.arange(len(data)) != i]

        # calculate the interesting quantity with the slightly reduced dataset.
        resampled_quantity[i] = interesting_quantity_func(resample)

    # find the error on the interesting quantity using the calculated values.
    error = np.sqrt(np.sum((resampled_quantity - all_samples) ** 2))
    return error


def get_bootstrap_error(data, interesting_quantity_func):
        """working bootstrap method."""

        # how many resamples to do. 1000 should work well.
        k = 1000

        # prepare array for quantities calculated with resampling.
        resampled_quantity = np.zeros(k)

        # resample k times.
        for i in range(k):

            # take a sample from the data, same length as the data set but
            # resampled with replacement.
            resample = np.random.choice(data, len(data))

            # calculate the interesting quantity with resampled dataset.
            resampled_quantity[i] = interesting_quantity_func(resample)

        # find the standard deviation of the newly calculated values.
        error = np.sqrt(np.var(resampled_quantity))
        return error


# get_all_data()
analyse()