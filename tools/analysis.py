from glob import glob
import matplotlib.pyplot as plt
import numpy as np
import os

class SimPlotter(object):
    def __init__(self, fname, num_cols, step = 1):
        self.fname = fname
        self.num_cols = num_cols
        self.step = step

        data = np.loadtxt(self.fname)
        self.times = data[::step, 0]
        self.populations = data[::step, 1::4]
        self.polarizations = data[::step, 3::4] + 1j*data[::step, 4::4]

    def population_plot(self):
        fig, axes = plt.subplots()

        for i in range(self.num_cols):
            axes.plot(self.times, self.populations[:,i])

        axes.set_xlabel("Time (ps)")
        axes.set_ylabel("Population")
        axes.set_ylim((-0.08,1.08))
        axes.set_title(self.fname)

        return fig



base_path = "/home/connor/Scratch/extensive"
data_paths = glob(os.path.join(base_path,"*_pi_pulse/sim/density*"))

sp = SimPlotter("/home/connor/Scratch/extensive/0.5_pi_pulse/sim/density_01/output.dat", 1024, 3)

fig = sp.population_plot()

plt.show()
