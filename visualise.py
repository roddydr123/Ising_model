import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from threading import Thread
from physics import Ising

class Animation(object):

    def __init__(self, system_size, temperature, method):
        # Set up the model
        self.model = Ising(system_size, temperature, method)
        self.nsteps = 50
        self.printFreq = 1

        # Set up matplotlib variables
        self.fig, self.ax = plt.subplots()
        
        # Draw the spins as an image in matplotlib
        self.implot = self.ax.imshow(self.model.spins)
        self.ani = None # For storing the animation object
        
    def run(self):
        # Create a separate thread to run the simulation model
        # Note that variables are shared between the threads
        thread = Thread(target=self.model.run,
                        args=(self.nsteps, self.printFreq,))

        # Start the simulation thread
        thread.start()

        # Start the animation in the main thread
        self.ani = animation.FuncAnimation(self.fig, self.animate,
                                           interval=10, blit=True)
        plt.show()

        # Stop the simulation after the animation window has closed
        self.model.stop()

    def animate(self, frame):
        # This function is called every time when the animation object updates
        # the screen. The argument frame (required) gives the index of the 
        # current animation frame
        
        # Update the image to show the latest configuration of the spins
        self.implot.set_data(self.model.spins)

        # Return the "artists" whose contents have been updated. Here we have 
        # changed the image data, so we need to return the artist that is 
        # responsible for drawing the image, which is implot
        return self.implot,


# Main entry point of the program
if __name__ == "__main__":
    
    # Read input arguments
    args = sys.argv
    if (len(args) != 4):
        print("Usage visualise.py system_size, temperature, method")
        sys.exit(1)

    system_size = tuple(args[1])
    temperature = float(args[2])
    method = str(args[3])

    # Set up and run the visualisation
    view = Animation(system_size, temperature, method)
    view.run()