import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import itertools as it
import sys


def updateSpinsGlauber(grid, rands, count, temperature, probs, J, grid_size):

    # choose a random site
    spin_indices = rands[count][:2]

    # calculate change in energy if flipped
    deltaE = get_deltaE(spin_indices, grid, J, grid_size)

    prob_flip = np.exp(-deltaE/temperature)

    if probs[count] < prob_flip:
        grid[spin_indices[0], spin_indices[1]] *= -1
    
    return grid



def updateSpinsKawasaki(grid, rands, count, temperature, probs, J, grid_size):

    # choose a random site
    first_spin_indices = rands[count][:2]
    second_spin_indices = rands[count][2:]

    first_spin_value = grid[first_spin_indices[0], first_spin_indices[1]]
    second_spin_value = grid[second_spin_indices[0], second_spin_indices[1]]

    if first_spin_value == second_spin_value:
        # if the spins are the same, don't swap them.
        return grid

    # calculate change in energy if flipped.
    deltaE = get_deltaE(first_spin_indices, grid, J, grid_size) + get_deltaE(second_spin_indices, grid, J, grid_size)

    # check if nearest neighbours by computing distance between the spins.
    if np.linalg.norm([first_spin_indices[0] - second_spin_indices[0], first_spin_indices[1], second_spin_indices[1]])%grid_size == 1:
        deltaE +=4

    prob_flip = np.exp(-deltaE/temperature)

    if probs[count] < prob_flip:
        grid[first_spin_indices[0], first_spin_indices[1]] *= -1
        grid[second_spin_indices[0], second_spin_indices[1]] *= -1

    return grid



def get_nn_spins(indices, grid, grid_size):
    """Find spins of nearest neighbours to spin at [indices].

    Args:
        indices (list): indices of the spin to be checked.

    Returns:
        array: The four nearest neighbour spins.
    """

    nn_spins = [grid[(indices[0] + 1)%grid_size, indices[1]],
                            grid[indices[0] - 1, indices[1]],
                            grid[indices[0], indices[1] - 1],
                            grid[indices[0], (indices[1] + 1)%grid_size]]
    return nn_spins


def get_deltaE(indices, grid, J, grid_size):
    """Finds the energy difference if a spin is flipped.

    Args:
        indices (list): indices of the spin to be (maybe) flipped.

    Returns:
        float: difference in energy between original state and spin flip.
    """

    spin_value = grid[indices[0], indices[1]]
    nn_spins = get_nn_spins(indices, grid, grid_size)

    return 2 * J * spin_value * np.sum(nn_spins)



def get_susceptibility(magnetisations=None):
    """Calculate magnetic susceptibility from all measurements.

    Returns:
        float: Magnetic susceptibility.
    """

    if magnetisations is None:
        magnetisations = magnetisation_list

    sus = 1/(temperature * grid_size**2) * np.var(magnetisations)
    return sus



def get_total_energy(grid, J):
    """Calculate the total energy of the system by summing over spins.

    Returns:
        float: Total energy
    """

    # # create a list of the spin indices, excluding ones which will be double counted.
    # indices = list(it.filterfalse(lambda x: x[0] >= x[1], np.ndindex(grid)))

    # # find the nearest neighbour spins.
    # nn_spins = list(map(lambda index: get_nn_spins((index[0], index[1])), indices))

    # # only take some of the spin matrix to avoid double counting.
    # new_grid = list(map(lambda index: grid[index[0], index[1]], indices))

    # # calculate total energy.
    # total_energy = -1 * J * np.sum(new_grid * np.sum(nn_spins, axis=1))
    total_energy = 0

    return total_energy



def get_heat_capacity(energies=None):
    """Calculate the heat capacity per spin.

    Args:
        energies (array, optional): sample of total energies. Defaults to None.
        no_states (int, optional): number of states to sample. Defaults to None.

    Returns:
        float: heat capacity
    """
    if energies is None:
        energies = energy_list

    C = 1/(temperature**2 * grid_size**2) * np.var(energies)
    return C



def get_bootstrap_error(n, k, q):
    """Compute error on heat capacity at the end of the run.

    Args:
        n (int): number of states to sample per sp heat calculation.
        k (int): number of times to repeat resampling.

    Returns:
        float: error on specific heat for this temperature.
    """

    quantity = []

    if q == "c":
        for i in range(k):
            energies = np.random.choice(energy_list, n)
            quantity.append(get_heat_capacity(energies=energies))
    elif q == "x":
        for i in range(k):
            mags = np.random.choice(magnetisation_list, n)
            quantity.append(get_susceptibility(magnetisations=mags))

    g = np.average(quantity)**2
    h = np.average(np.array(quantity)**2)
    error = np.sqrt(h-g)
    return error


def run_sim(kT, update_func, nsweeps, grid_size, vis, grid, J):

    # if we wanted a visualisation, make a figure for it.
    if vis:
        fig, ax = plt.subplots()
        im = ax.imshow(grid, animated=True, cmap='bwr')

    # enough for 1 temperature: sweeps * iterations per sweep
    n_randints = nsweeps * grid_size**2
    # random indices to choose spins
    rands = np.random.randint(0, grid_size, size=(n_randints, 4))
    # random numbers to compare probabilities with
    probs = np.random.rand(n_randints)

    wait = 100

    for n in tqdm(range(nsweeps)):
        for i in range(grid_size**2):
            grid = update_func(grid, rands, (n * grid_size**2 + i), kT, probs, J, grid_size)
            # print(grid)

        # every 50 sweeps update the animation.
        if n % 10 == 0 and vis:
            
            plt.cla()
            im = ax.imshow(grid, animated=True, cmap='bwr', interpolation='bicubic')
            plt.draw()
            plt.pause(0.00001)

        # every 10 sweeps record the free energy.
        if n % 10 == 0:

            magnetisation = np.sum(grid)
            total_energy = get_total_energy(grid, J)

            with open("ising.dat", "a") as f:
                f.write(f"{n}, {magnetisation}, {total_energy}\n")

        # check if the simulation has converged.
        if n % 10 == 0 and n > wait and update_func == updateSpinsGlauber:
            data = np.genfromtxt('ising.dat', delimiter=',', skip_header=1)
            if len(set(np.round(data[-10:-1][:,1], 7))) == 1:
                # convergence
                break


def main():

    try:
        _, vis, grid_size, temperature, method = sys.argv
    except:
        print("Usage ising_model.py <vis> <grid_size> <temperature> <method>")
        sys.exit()


    grid_size = int(grid_size)
    temperature = float(temperature)

    if method == "G":
        update_func = updateSpinsGlauber
    elif method == "K":
        update_func = updateSpinsKawasaki
    else:
        raise TypeError("Check spins matrix type")
    

    # random grid
    grid = np.random.choice([-1, 1], [grid_size, grid_size])
    
    # show a visualisation or not.
    if vis == "vis":
        vis = True
    elif vis == "novis":
        vis = False
    else:
        sys.exit()

    J = 1.0
    k = 1.0
    nsweeps = 1000

    # clear a file to write free energy data to.
    with open("ising.dat", "w") as f:
        f.write("iteration, magnetisation, total energy\n")

    run_sim(k * temperature, update_func, nsweeps, grid_size, vis, grid, J)


main()