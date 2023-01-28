import os
import sys
import numpy as np


J = 1


class Ising(object):

    def __init__(self, system_size, temperature, method):
        self.system_size = system_size
        self.spins = np.zeros((system_size, system_size), dtype=float)

        # make spins random
        self.spins = np.random.rand(system_size, system_size)
        self.spins[self.spins > 0.5] = 1
        self.spins[self.spins < 0.5] = -1

        self.temperature = temperature
        self.method = method

        self.stopSim = False # Flag to stop the simulation thread

    def updateSpinsGlauber(self):
        # choose a random site
        spin_index = (np.random.randint(self.system_size), np.random.randint(self.system_size))

        # calculate change in energy if flipped
        deltaE = self.get_deltaE(spin_index)

        prob_flip = np.exp(-deltaE/self.temperature)

        prob_rand = np.random.random()

        if prob_rand < prob_flip:
            self.spins[spin_index[0], spin_index[1]] *= -1
        return

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

        prob_flip = np.exp(-deltaE/self.temperature)

        prob_rand = np.random.random()

        if prob_rand < prob_flip:
            self.spins[first_spin_indices[0], first_spin_indices[1]] *= -1
            self.spins[second_spin_indices[0], second_spin_indices[1]] *= -1
        return

    def printSpins(self):

        outfile = "spins.dat"
        # Write the lattice first to a temporary file
        with open(outfile+".tmp", "w") as writer:
            for i in range(self.system_size):
                for j in range(self.system_size):
                    writer.write(f"{i} {j} {self.spins[i,j]}\n")
                writer.write("\n")
        # Rename the temporary file to the output file
        os.rename(outfile+".tmp", outfile)

    def run(self, nsteps, printFreq):
        self.stopSim = False

        for n in range(nsteps):
            for M in range(self.system_size):
                for N in range(10):
                    if (self.stopSim): break

                    if self.method == "G":
                        self.updateSpinsGlauber()
                    elif self.method == "K":
                        self.updateSpinsKawasaki()
                    else:
                        print("Method: G for Glauber or K for Kawasaki")
                        sys.exit(1)

            if (n % printFreq == 0):
                # print(n)
                self.printSpins()

    def stop(self):
        self.stopSim = True

    def get_total_magnetisation(self):
        return np.sum(self.spins)

    def get_total_energy(self):
        return 0


# Main entry point of the program without visualisation
if __name__ == "__main__":

    # Read input arguments
    args = sys.argv
    if (len(args) != 4):
        print("Usage visualise.py system_size, temperature, method")
        sys.exit(1)

    system_size = int(args[1])
    temperature = float(args[2])
    method = str(args[3])
    nsteps = 10
    printFreq = 1
    
    # Set up and run the model
    model = Ising(system_size, temperature, method)
    model.run(nsteps, printFreq)
