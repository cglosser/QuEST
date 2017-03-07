import matplotlib.pyplot as plt
import numpy as np
from itertools import chain, zip_longest

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)

lo_alpha, hi_alpha = 0.075, 1

path = "/home/connor/Scratch/neighborhood/full_neighborhood/sim/neighborhood_03/output.dat"

output = np.loadtxt(path)

times = output[::3,0]
pop = output[::3,1::4]

lines = []
for i in range(1024):
    lines.append(plt.plot(times, pop[:,i], alpha=lo_alpha, zorder=1024 - i))

plt.show(block=False)


interesting_idx = set()

for group in grouper(enumerate(chain.from_iterable(lines)), 10):
    for idx, line in group:
        line.set_alpha(hi_alpha)

    plt.draw()

    query = input("Any interesting?")
    if query.lower() in ['y', 'yes']:
        for idx, line in group:
            line.set_alpha(lo_alpha)

        for idx, line in group:
            line.set_alpha(hi_alpha)
            plt.draw()

            query = input("Interesting ({})? ".format(idx))
            if query.lower() in ['y', 'yes']:
                interesting_idx.add(idx)
                line.set_alpha(lo_alpha)
            else:
                line.remove()
    else:
        for idx, line in group:
            line.remove()


with open("interesting.dat", 'w') as f:
    for idx in interesting_idx:
        f.write(str(idx) + "\n")
