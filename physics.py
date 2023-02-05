import numpy as np
import itertools as it


J = 1


class Ising(object):



    def __init__(self, system_size, temperature, method, nstep, spins=None):
        np.random.seed(5)
        self.system_size = system_size

        # make system of random spins up or down if kawasaki or chosen
        if (spins == "random" or method == "K"):
            self.spins = np.random.choice([-1, 1], [system_size, system_size])
        # if glauber check if spins passed in.
        elif (spins is not None and method == "G"):
            self.spins = spins
        # if glauber set spins all to zero at start.
        else:
            self.spins = np.ones((system_size, system_size))

        self.temperature = temperature
        self.method = method
        self.magnetisation_list = []
        self.energy_list = []

        # enough for 1 temperature: sweeps * iterations per sweep
        n_randints = nstep * system_size**2

        # random indices to choose spins
        self.rands = np.random.randint(0, system_size, size=(n_randints, 4))

        # random numbers to compare probabilities with
        self.probs = np.random.rand(n_randints)

        # keeps track of which iteration we're on for indexing random arrays.
        self.count = 0


    def updateSpinsGlauber(self):
        # choose a random site

        spin_indices = self.rands[self.count][:2]

        # calculate change in energy if flipped
        deltaE = self.get_deltaE(spin_indices)

        prob_flip = np.exp(-deltaE/self.temperature)

        if self.probs[self.count] < prob_flip:
            self.spins[spin_indices[0], spin_indices[1]] *= -1



    def updateSpinsKawasaki(self):
        # choose a random site
        first_spin_indices = self.rands[self.count][:2]
        second_spin_indices = self.rands[self.count][2:]

        first_spin_value = self.spins[first_spin_indices[0], first_spin_indices[1]]
        second_spin_value = self.spins[second_spin_indices[0], second_spin_indices[1]]

        if first_spin_value == second_spin_value:
            # if the spins are the same, don't swap them.
            return

        # calculate change in energy if flipped.
        deltaE = self.get_deltaE(first_spin_indices) + self.get_deltaE(second_spin_indices)

        # check if nearest neighbours by computing distance between the spins.
        if np.linalg.norm([first_spin_indices[0] - second_spin_indices[0], first_spin_indices[1], second_spin_indices[1]])%self.system_size == 1:
            deltaE +=4

        prob_flip = np.exp(-deltaE/self.temperature)

        # prob_rand = np.random.random()

        if self.probs[self.count] < prob_flip:
            self.spins[first_spin_indices[0], first_spin_indices[1]] *= -1
            self.spins[second_spin_indices[0], second_spin_indices[1]] *= -1



    def get_nn_spins(self, indices):
        """Find spins of nearest neighbours to spin at [indices].

        Args:
            indices (list): indices of the spin to be checked.

        Returns:
            array: The four nearest neighbour spins.
        """

        nn_spins = [self.spins[(indices[0] + 1)%self.system_size, indices[1]],
                             self.spins[indices[0] - 1, indices[1]],
                             self.spins[indices[0], indices[1] - 1],
                             self.spins[indices[0], (indices[1] + 1)%self.system_size]]
        return nn_spins


    def get_deltaE(self, indices):
        """Finds the energy difference if a spin is flipped.

        Args:
            indices (list): indices of the spin to be (maybe) flipped.

        Returns:
            float: difference in energy between original state and spin flip.
        """

        spin_value = self.spins[indices[0], indices[1]]
        nn_spins = self.get_nn_spins(indices)
        return 2 * J * spin_value * np.sum(nn_spins)



    def get_total_magnetisation(self):
        magnetisation = np.sum(self.spins)
        self.magnetisation_list.append(magnetisation)
        return magnetisation



    def get_susceptibility(self, magnetisations=None):
        """Calculate magnetic susceptibility from all measurements.

        Returns:
            float: Magnetic susceptibility.
        """

        if magnetisations is None:
            magnetisations = self.magnetisation_list

        sus = 1/(self.temperature * self.system_size**2) * np.var(magnetisations)
        return sus



    def get_total_energy(self):
        """Calculate the total energy of the system by summing over spins.

        Returns:
            float: Total energy
        """

        spins = self.spins

	    # create a list of the spin indices, excluding ones which will be double counted.
        indices = list(it.filterfalse(lambda x: x[0] >= x[1], np.ndindex(spins.shape)))

	    # find the nearest neighbour spins.
        nn_spins = list(map(lambda index: self.get_nn_spins((index[0], index[1])), indices))

	    # only take some of the spin matrix to avoid double counting.
        new_spins = list(map(lambda index: spins[index[0], index[1]], indices))

	    # calculate total energy.
        total_energy = -1 * J * np.sum(new_spins * np.sum(nn_spins, axis=1))

        self.energy_list.append(total_energy)
        return total_energy



    def get_heat_capacity(self, energies=None):
        """Calculate the heat capacity per spin.

        Args:
            energies (array, optional): sample of total energies. Defaults to None.
            no_states (int, optional): number of states to sample. Defaults to None.

        Returns:
            float: heat capacity
        """
        if energies is None:
            energies = self.energy_list

        C = 1/(self.temperature**2 * self.system_size**2) * np.var(energies)
        return C



    def get_bootstrap_error(self, n, k, q):
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
                energies = np.random.choice(self.energy_list, n)
                quantity.append(self.get_heat_capacity(energies=energies))
        elif q == "x":
            for i in range(k):
                mags = np.random.choice(self.magnetisation_list, n)
                quantity.append(self.get_susceptibility(magnetisations=mags))

        g = np.average(quantity)**2
        h = np.average(np.array(quantity)**2)
        error = np.sqrt(h-g)
        return error
