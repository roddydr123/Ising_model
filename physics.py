import numpy as np


J = 1


class Ising(object):

    def __init__(self, system_size, temperature, method):
        self.system_size = system_size

        # make system of random spins up or down
        self.spins = np.random.choice([-1, 1], [system_size, system_size])

        self.temperature = temperature
        self.method = method
        self.magnetisation_list = []

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

    def get_deltaE(self, indices):

        spin_value = self.spins[indices[0], indices[1]]

        added = list(indices)
        if indices[0] == self.system_size - 1:
            added[0] = -1
        if indices[1] == self.system_size - 1:
            added[1] = -1

        nn_spins = np.array([self.spins[added[0] + 1, indices[1]], self.spins[indices[0] - 1, indices[1]],
                    self.spins[indices[0], indices[1] - 1], self.spins[indices[0], added[1] + 1]])

        return 2 * J * spin_value * np.sum(nn_spins)

    def get_total_magnetisation(self):
        magnetisation = np.sum(self.spins)
        self.magnetisation_list.append(magnetisation)
        return magnetisation

    def get_susceptibility(self):
        # (average of M)**2
        g = np.average(self.magnetisation_list)**2
        # average of (M**2)
        h = np.average(np.array(self.magnetisation_list)**2)
        n = h - g
        return n/((self.system_size)**2 * self.temperature)

    def get_total_energy(self):
        return 0
