#ifndef HISTORY_H
#define HISTORY_H

#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 3>;
}

template <class soltype>
class Integrator::History {
 public:
  History(const int, const int, const int);
  soltype_array<soltype> array;

 private:
};

template <class soltype>
Integrator::History<soltype>::History(const int num_particles, const int window,
                                      const int num_timesteps)
    : array(boost::extents[num_particles][
          typename soltype_array<soltype>::extent_range(-window, num_timesteps)]
                            [2])

{
}

#endif
