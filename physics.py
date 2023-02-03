import numpy as np


J = 1


class Ising(object):



    def __init__(self, system_size, temperature, method):
        np.random.seed(10)
        self.system_size = system_size

        # make system of random spins up or down
        self.spins = np.random.choice([-1, 1], [system_size, system_size])

        self.temperature = temperature
        self.method = method
        self.magnetisation_list = []
        self.energy_list = []

        # # random indices to choose spins
        # self.rands = np.random.randint(0, system_size, size=(n_randints, 4))

        # # random numbers to compare probabilities with
        # self.probs = 3



    def updateSpinsGlauber(self):
        # choose a random site
        spin_indices = (np.random.randint(self.system_size), np.random.randint(self.system_size))

        # calculate change in energy if flipped
        deltaE = self.get_deltaE(spin_indices)

        prob_flip = np.exp(-deltaE/self.temperature)

        prob_rand = np.random.random()

        if prob_rand < prob_flip:
            self.spins[spin_indices[0], spin_indices[1]] *= -1



    def updateSpinsKawasaki(self):
        # choose a random site
        first_spin_indices = (np.random.randint(self.system_size), np.random.randint(self.system_size))
        second_spin_indices = (np.random.randint(self.system_size), np.random.randint(self.system_size))

        first_spin_value = self.spins[first_spin_indices[0], first_spin_indices[1]]
        second_spin_value = self.spins[second_spin_indices[0], second_spin_indices[1]]

        if first_spin_value == second_spin_value:
            # if the spins are the same, don't swap them.
            return

        # calculate change in energy if flipped
        deltaE = self.get_deltaE(first_spin_indices) + self.get_deltaE(second_spin_indices)

        if np.linalg.norm([first_spin_indices[0] - second_spin_indices[0], first_spin_indices[1], second_spin_indices[1]])%self.system_size == 1:
            deltaE +=4

        prob_flip = np.exp(-deltaE/self.temperature)

        prob_rand = np.random.random()

        if prob_rand < prob_flip:
            self.spins[first_spin_indices[0], first_spin_indices[1]] *= -1
            self.spins[second_spin_indices[0], second_spin_indices[1]] *= -1



    def get_nn_spins(self, indices):
        """Find spins of nearest neighbours to spin at [indices].

        Args:
            indices (list): indices of the spin to be checked.

        Returns:
            array: The four nearest neighbour spins.
        """

        nn_spins = np.array([self.spins[(indices[0] + 1)%self.system_size, indices[1]],
                             self.spins[indices[0] - 1, indices[1]],
                             self.spins[indices[0], indices[1] - 1],
                             self.spins[indices[0], (indices[1] + 1)%self.system_size]])
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



    def get_susceptibility(self, magnetisations=None, no_states=None):
        """Calculate magnetic susceptibility from all measurements.

        Returns:
            float: Magnetic susceptibility.
        """

        if magnetisations is None:
            magnetisations = self.magnetisation_list
        if no_states is None:
            no_states = self.system_size

        # (average of M)**2
        # g = np.average(magnetisations)**2
        # # average of (M**2)
        # h = np.average(np.array(magnetisations)**2)
        # n = h - g
        norm = 1/(self.system_size**2 * self.temperature)
        sus = norm * (np.average(np.square(magnetisations)) - np.square(np.average(magnetisations)))
        return sus
        # return n/((no_states)**2 * self.temperature)



    def get_total_energy(self):
        """Calculate the total energy of the system by summing over spins.

        Returns:
            float: Total energy
        """
        total_energy = 0

        for i in range(self.system_size):
            # don't repeat sum
            for j in range(i+1, self.system_size):
                nn_spins = self.get_nn_spins([i,j])
                total_energy += self.spins[i,j] * np.sum(nn_spins)

        total_energy *= -1 * J
        self.energy_list.append(total_energy)
        return total_energy



    def get_heat_capacity(self, energies=None, no_states=None):
        """Calculate the heat capacity per spin.

        Args:
            energies (array, optional): sample of total energies. Defaults to None.
            no_states (int, optional): number of states to sample. Defaults to None.

        Returns:
            float: heat capacity
        """
        if energies is None:
            energies = self.energy_list
        if no_states is None:
            no_states = self.system_size

        norm_fact = 1 / (self.system_size**2 * self.temperature**2)
        heat_cap = norm_fact * (np.average(np.square(energies)) - np.square(np.average(energies)))

        return heat_cap

        # g = np.average(energies)**2

        # h = np.average(np.array(energies)**2)
        # n = h-g
        # return np.var(energies) / (no_states**2 * self.temperature**2)
        # return n/((no_states)**2 * self.temperature**2)



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
                quantity.append(self.get_heat_capacity(energies=energies, no_states=n))
        elif q == "x":
            for i in range(k):
                mags = np.random.choice(self.magnetisation_list, n)
                quantity.append(self.get_susceptibility(magnetisations=mags, no_states=n))

        g = np.average(quantity)**2
        h = np.average(np.array(quantity)**2)
        error = np.sqrt(h-g)
        return error
