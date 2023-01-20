import os
import sys
import numpy as np


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

        self.dt = 0.1 # Timestep
        self.stopSim = False # Flag to stop the simulation thread

    
    def updateSpins(self):
        # choose a random site
        spin_index = (np.random.randint(self.system_size), np.random.randint(self.system_size))
        # spin_index = (3, 3)
        spin_value = self.spins[spin_index[0], spin_index[1]]

        # find nearest neighbour spins. Check if adding 1 will be outside the matrix.
        added = list(spin_index)
        if spin_index[0] == system_size - 1:
            added[0] = -1
        if spin_index[1] == system_size - 1:
            added[1] = -1

        nn_spins = np.array([self.spins[added[0] + 1, spin_index[1]], self.spins[spin_index[0] - 1, spin_index[1]],
                    self.spins[spin_index[0], spin_index[1] - 1], self.spins[spin_index[0], added[1] + 1]])
    
        # print(nn_spins)
        print(self.spins)
        # print(spin_value)

        # calculate change in energy if flipped

        deltaE = -1 * np.sum(spin_value * nn_spins)
        print(deltaE)
        print(spin_index)

        prob_flip = 1
        if deltaE > 0:
            prob_flip = np.exp(-deltaE/self.temperature)

        prob_rand = np.random.rand(1)

        if prob_rand > 1 - prob_flip:
            self.spins[spin_index[0], spin_index[1]] *= -1

        # just now they always get flipped...


        


    def printSpins(self):
        outfile = "spins.dat"
        # Write the lattice first to a temporary file
        with open(outfile+".tmp", "w") as writer:
            for i in range(self.system_size):
                for j in range(self.system_size):
                    writer.write("{:d} {:d} {:d}\n".format(
                        i, j, self.spins[i,j]))
                writer.write("\n")
        # Rename the temporary file to the output file
        os.rename(outfile+".tmp", outfile)

    
    def run(self, nsteps, printFreq):
        self.stopSim = False
        for i in range(nsteps):
            if (self.stopSim): break
            self.updateSpins()
            
            # Draw and output the lattice configuration every printFreq steps
            if (i % printFreq == 0):
                print("Step {:d}".format(i))
                # self.printSpins()


    def stop(self):
        self.stopSim = True


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
    nsteps = 5
    printFreq = 1
    
    # Set up and run the model
    model = Ising(system_size, temperature, method)
    model.run(nsteps, printFreq)
