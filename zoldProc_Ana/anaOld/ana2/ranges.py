import numpy as np
import matplotlib.pyplot as plt

print("calculating gaseous and liquid ranges")

plt.ion()

cEn, cRange = np.loadtxt("./dedx_carbon.txt", skiprows = 5, delimiter = "|", usecols = (0, 5), unpack = True)


# np.loadtxt("dedx_fluorine.txt", delimiter = "|")

fig, axes = plt.subplots(1, 1, figsize = (5, 5))

plt.show()