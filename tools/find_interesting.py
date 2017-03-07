import matplotlib.pyplot as plt
import numpy as np
from itertools import chain

lo_alpha, hi_alpha = 0.075, 1

path = "/home/connor/Scratch/neighborhood/full_neighborhood/sim/neighborhood_00/output.dat"

output = np.loadtxt(path)

times = output[::3,0]
pop = output[::3,1::4]

lines = []
for i in range(1024):
    lines.append(plt.plot(times, pop[:,i], alpha=lo_alpha))

plt.show(block=False)


interesting_idx = set()

for idx, line in enumerate(chain.from_iterable(lines)):
    line.set_alpha(hi_alpha)
    plt.draw()

    query = input("Interesting ({})? ".format(idx))
    if query.lower() in ['y', 'yes']:
        interesting_idx.add(idx)
        line.set_alpha(lo_alpha)
    else:
        line.remove()


with open("interesting.dat", 'w') as f:
    for idx in interesting_idx:
        f.write(str(idx) + "\n")
