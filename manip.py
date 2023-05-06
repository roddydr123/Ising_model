import numpy as np
from physics import Ising


n = 200
k = 50

temps = np.round(np.arange(1, 3, 0.1), 3)[:1]

for temp in temps:

    model = Ising(50, temp, "G", 10100)
    data = np.genfromtxt(f"data/{temp}.mags.G4.dat")

    # energies_list = data
    magnetisations_list = data

    quantity = []

    # for i in range(k):
    #     mags = np.random.choice(magnetisation_list, n)
    #     quantity.append(model.get_susceptibility(magnetisations=mags))

    # g = np.average(quantity)**2
    # h = np.average(np.array(quantity)**2)
    # error = np.sqrt(h-g)

    print(np.std(magnetisations_list))
    print(np.std(abs(magnetisations_list)))

    print(data)






# av_energy = np.average(energies)
# av_mag = np.average(mags)

# print(av_energy, av_mag)

# norm = 1/(50**2 * 1)
# sus = norm * (np.average(np.square(mags)) - np.square(np.average(mags)))
# C = norm * (np.average(np.square(energies)) - np.square(np.average(energies)))

# C = 1/(temperature * 2500) * (np.var(energies))

# sus = 1/(temperature * 2500) * (np.var(mags))

# print(sus, C)