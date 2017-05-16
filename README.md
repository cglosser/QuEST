# Quantum Electromagnetics Simulation Toolbox (QuEST)

Simulation software to model the evolution & electromagnetic interactions of two-level quantum dots. Designed and built at Michigan State University.

## Prerequisites

* C++14-compatible compiler (tested with [gcc](https://gcc.gnu.org/) and [clang](https://clang.llvm.org/))
* [Eigen3](http://eigen.tuxfamily.org) (at least v3.2.2)
* Boost (at least v.1.55.0):
  * [program options](http://www.boost.org/doc/libs/1_55_0/doc/html/program_options.html)
  * [multiarray](http://www.boost.org/doc/libs/1_55_0/libs/multi_array/doc/index.html)
  * [test](http://www.boost.org/doc/libs/1_64_0/libs/test/doc/html/index.html)
* [SILO](https://wci.llnl.gov/simulation/computer-codes/silo) (optional)

## Building

QuEST relies on [CMake](https://cmake.org/) to generate appropriate compile scripts. To build the executable, first run

    mkdir build
    cd build

followed by

    cmake $PATH_TO_QUEST && make

This will attempt to build

* `quest` (simulation executable)
* `qtest` (unit test executable)
* `point_gen` (utility to quickly generate distributions of points)
* `siloify` (utility to convert `quest` output to the [SILO](https://wci.llnl.gov/simulation/computer-codes/silo) file format for use in e.g. [VisIt](https://wci.llnl.gov/simulation/computer-codes/visit/)

## Running

QuEST requires three files to specify the parameters of a simulation:

1. `input.cfg`
2. `dots.cfg`
3. `pulse.cfg`

Of these, `input.cfg` has the most flexibility; you can specify alternate paths
to the other configuration files as well as modify the simulation parameters
(speed of light, number of particles, timestep, etc.). Run `./quest --help` for
details of the possible parameter options.

`dots.cfg` contains a list of quantum dots, one-per-line, each with the
following format:

    x y z omega_0 T1 T2 dx dy dz

`pulse.cfg` specifies the incident pulse(s) with the following format:

    E0 delay sigma omega_L kx ky kz px py pz


