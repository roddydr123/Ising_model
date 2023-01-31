from physics import Ising
import numpy as np
from tqdm import tqdm

J=1.0
nstep=10100
wait=100

lx=50
ly=lx

method = "G"

def onerun(kT):
    # no visualisation
    # run a sim for 10000 sweeps
    # start measurements after 100 sweeps
    # record the magnetisation and susceptibility every 10 sweeps

    model = Ising(lx, kT, method)

    runfile = f"data/{kT}output.dat"

    # clear file
    f = open(runfile, "w")
    f.close()

    for n in tqdm(range(nstep)):
        for i in range(lx):
            for j in range(ly):
                if method == "G":
                    model.updateSpinsGlauber()
                if method == "K":
                    model.updateSpinsKawasaki()

        if n % 10 == 0 and n > wait:

            total_magnetisation = model.get_total_magnetisation()
            total_energy = model.get_total_energy()
            susceptibility = model.get_susceptibility()
            heat_capacity = model.get_heat_capacity()
            
            # record required data
            with open(runfile, "a") as file:
                file.write(f"{total_magnetisation}, {total_energy}, {susceptibility}, {heat_capacity}\n")

    # compute average specific heat for this temp
    # read totals from file
    read_data = np.genfromtxt(runfile, delimiter=",")

    av_sp_heat = np.average(read_data[:,3])
    av_energy = np.average(read_data[:,1])
    av_magnetisation = np.average(read_data[:,0])
    av_susceptibility = np.average(read_data[:,2])

    c_error = model.get_bootstrap_error(200, 50)
    with open("data/bs_errors.dat", "a") as errorfile:
        errorfile.write(f"{kT}, {av_sp_heat}, {c_error}, {av_energy}, {av_magnetisation}, {av_susceptibility}\n")


def main():

    # clear file
    f = open("data/bs_errors.dat", "w")
    f.close()

    # make array of all the temperatures to loop through
    temps = np.round(np.arange(1, 3, 0.1), 3)

    for temp in tqdm(temps):
        onerun(temp)


main()