#ifndef UNIVERSE_H
#define UNIVERSE_H

#include <boost/program_options.hpp>

struct Universal{
  size_t num_particles;
  double c0, hbar;
  double simulation_time;
};

extern Universal Universe;

#endif
