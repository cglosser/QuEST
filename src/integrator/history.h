#ifndef HISTORY_H
#define HISTORY_H

#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

#include "../magnetic_particle.h"

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 3>;
}

template <class soltype>
class Integrator::History {
 public:
  static constexpr int DERIV_0 = 0, DERIV_1 = 1;

  History(const int, const int, const int);
  soltype_array<soltype> array;

  void fill(const soltype &);
  void initialize_past(const std::shared_ptr<const DotVector> &);

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

template <class soltype>
void Integrator::History<soltype>::fill(const soltype &val)
{
  std::fill(array.data(), array.data() + array.num_elements(), val);
}

template <class soltype>
void Integrator::History<soltype>::initialize_past(const std::shared_ptr<
                                                   const DotVector> &dots)
{
  for(int time_idx=array.index_bases()[1]; time_idx <= 0; ++time_idx) {
    for(int i=0; i < static_cast<int>(array.shape()[0]); ++i)
      array[i][time_idx][DERIV_0] = (*dots)[i].magnetization();
  }
}

#endif
