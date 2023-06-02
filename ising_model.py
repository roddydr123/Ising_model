import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
import sys
from scipy.signal import convolve2d


def updateSpinsGlauber(grid, rands, count, temperature, probs, J, grid_size):
    """
    One step of the Glauber algorithm. Choose a spin, find the change in energy if flipped,
    definitely flip the spin if it would lower the energy of the system, or flip with
    Boltzmann probability if it would increase the energy of the system.

    Modifies the grid in place.
    """

    # choose a random site
    spin_indices = rands[count][:2]

    # calculate change in energy if flipped
    deltaE = get_deltaE(spin_indices, grid, J, grid_size)

    prob_flip = np.exp(-deltaE/temperature)

    if probs[count] < prob_flip:
        grid[spin_indices[0], spin_indices[1]] *= -1
    
    return grid



def updateSpinsKawasaki(grid, rands, count, temperature, probs, J, grid_size):
    """
    One step of the Kawasaki algorithm. Choose two spins randomly, find the change in
    energy if they were swapped, definitely swap them if it would lower the energy of
    the system, or swap with Boltzmann probability if it would increase the energy of
    the system. An extra correction is added if the spins are nearest neighbours.

    Modifies the grid in place.
    """

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

    # check if nearest neighbours by computing distance between the spins. Add a correction if they are nns.
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
    """Finds the energy difference if a spin is flipped. Muliply by 2 as flipping
    a spin changes the relationship by two J (e.g. if we used to have a -J and now
    have a +J, this is a difference of 2J).

    Args:
        indices (list): indices of the spin to be (maybe) flipped.

    Returns:
        float: difference in energy between original state and spin flip.
    """

    spin_value = grid[indices[0], indices[1]]
    nn_spins = get_nn_spins(indices, grid, grid_size)

    return 2 * J * spin_value * np.sum(nn_spins)


def get_total_energy(grid, J):
    """Calculate the total energy of the system by summing over spins.

    Returns:
        float: Total energy
    """

    # define the convolution kernel
    kernel = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]])

    # compute the sum of nearest neighbors with periodic boundary conditions
    neighbour_sum = convolve2d(grid, kernel, mode='same', boundary='wrap')

    total_energy = -J * (1/2) * np.sum(neighbour_sum * grid)

    return total_energy


def run_sim(kT, update_func, nsweeps, grid_size, vis, grid, J, filename):

    # if we wanted a visualisation, make a figure for it.
    if vis:
        fig, ax = plt.subplots()
        im = ax.imshow(grid, animated=True, cmap='bwr')


    # clear a file to write data to.
    with open(filename, "w") as f:
        f.write("iteration, magnetisation, total energy\n")

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

            with open(filename, "a") as f:
                f.write(f"{n}, {magnetisation}, {total_energy}\n")

        # check if the simulation has converged.
        if n % 10 == 0 and n > wait and update_func == updateSpinsGlauber:
            data = np.genfromtxt(filename, delimiter=',', skip_header=1)
            if len(set(np.round(data[-10:-1][:,1], 7))) == 1:
                # convergence
                break

    return grid


def main():

    try:
        _, vis, grid_size, temperature, method = sys.argv
    except:
        print("Usage ising_model.py <vis> <grid_size> <temperature> <method>")
        sys.exit()

    filename = input("Write to file: <name>|'d': ")
    if filename == "d":
        filename = "ising.dat"
    else:
        filename = f"{filename}.dat"

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

    run_sim(k * temperature, update_func, nsweeps, grid_size, vis, grid, J, filename)


if __name__=="__main__":
    main()