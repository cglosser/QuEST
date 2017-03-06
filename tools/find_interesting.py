import matplotlib.pyplot as plt
import numpy as np
from itertools import chain

lo_alpha, hi_alpha = 0.05, 1

path = "/home/connor/Scratch/neighborhood/full_neighborhood/sim/neighborhood_00/output.dat"

output = np.loadtxt(path)

times = output[:,0]
pop = output[:,1::3]

lines = []
for i in range(1024):
    lines.append(plt.plot(times, pop[:,i], alpha=lo_alpha))

plt.show(block=False)

for line in chain.from_iterable(lines):
    line.set_alpha(hi_alpha)
    plt.draw()

    input("Interesting? ")

    line.set_alpha(lo_alpha)
