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

    def polarization_plot(self):
        fig, (axes_abs, axes_arg) = plt.subplots(2, sharex = True)

        for i in range(self.num_cols):
            axes_abs.plot(self.times, np.abs(self.polarizations[:,i]))
            axes_arg.plot(self.times, np.angle(self.polarizations[:,i]))

        axes_abs.set_ylabel("Magnitude")
        axes_abs.set_title(self.fname)
        axes_abs.set_ylim((-0.08,0.54))

        axes_arg.set_xlabel("Time (ps)")
        axes_arg.set_ylabel("Phase")

        return fig

    def edwards_anderson_plot(self):
        fig, axes = plt.subplots()

        dthetas = np.angle(self.polarizations) - np.angle(self.polarizations[0])
        axes.plot(self.times, np.abs(np.mean(np.exp(1j*dthetas), axis = 1)))

        axes.set_xlabel("Time (ps)")
        axes.set_ylabel("Edwards-Anderson order parameter")
        axes.set_ylim((-0.08,1.08))
        axes.set_title(self.fname)

        return fig

    def phase_order_plot(self):
        fig, axes = plt.subplots()

        sinusoids = np.mean(np.exp(1j*np.angle(self.polarizations)), axis = 1)
        axes.plot(self.times, np.real(sinusoids), label = "real part")
        axes.plot(self.times, np.imag(sinusoids), label = "imag part")

        axes.set_xlabel("Time (ps)")
        axes.set_ylabel("Phase order parameter")
        axes.set_ylim((-1.08,1.08))
        axes.set_title(self.fname)
        axes.legend(loc="upper right")

        return fig

    def inverse_participation_plot(self):
        fig, axes = plt.subplots()

        axes.plot(self.times, np.sum(np.abs(self.polarizations)**4, axis = 1)/np.sum(np.abs(self.polarizations)**2, axis = 1)**2)

        axes.set_xlabel("Time (ps)")
        axes.set_ylabel("Inverse participation ratio")
        axes.set_title(self.fname)

        return fig

    def save_plots(self, path):
        self.population_plot().savefig(os.path.join(path, "population.png"))
        self.polarization_plot().savefig(os.path.join(path, "polarization.png"))
        self.edwards_anderson_plot().savefig(os.path.join(path, "edwards-anderson.png"))
        self.phase_order_plot().savefig(os.path.join(path, "phase_order.png"))
        self.inverse_participation_plot().savefig(os.path.join(path, "inverse_participation_ratio.png"))

def main():
    base_path = "/home/connor/Scratch/extensive"
    data_paths = sorted(glob(os.path.join(base_path,"*_pi_pulse/sim/density*")))

    for source in data_paths:
        print(source)

        output_dir = os.path.join(source, "figures")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        sp = SimPlotter(os.path.join(source, "output.dat"), 1024, 3)
        sp.save_plots(output_dir)

        plt.close('all')

if __name__=="__main__":
    main()
