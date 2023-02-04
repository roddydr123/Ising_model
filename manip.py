import numpy as np

temperature = 1


data = np.genfromtxt("outfile2.4.dat")

energies = data.T[1]
mags = data.T[0]

av_energy = np.average(energies)
av_mag = np.average(mags)

print(av_energy, av_mag)

norm = 1/(50**2 * 1)
sus = norm * (np.average(np.square(mags)) - np.square(np.average(mags)))
C = norm * (np.average(np.square(energies)) - np.square(np.average(energies)))

C = 1/(temperature * 2500) * (np.var(energies))

sus = 1/(temperature * 2500) * (np.var(mags))

print(sus, C)

H_errs.append(np.std( np.square(np.array(E_vals) - e_avg) / (no_particles * T**2)) / len(E_vals))