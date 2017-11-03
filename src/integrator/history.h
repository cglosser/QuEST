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

  inline namespace enums {
    enum DIMENSION { PARTICLES, TIMES, DERIVATIVES };
    enum ORDER { DERIV_0, DERIV_1 };
  }
}

template <class soltype>
class Integrator::History {
 public:
  History(const int, const int, const int, const int = 2);
  soltype_array<soltype> array;

  void fill(const soltype &);
  void initialize_past(const soltype &);

 private:
};

template <class soltype>
Integrator::History<soltype>::History(const int num_particles,
                                      const int window,
                                      const int num_timesteps,
                                      const int num_derivatives)
    : array(boost::extents[num_particles][
          typename soltype_array<soltype>::extent_range(-window, num_timesteps)]
                          [num_derivatives])

{
}

template <class soltype>
void Integrator::History<soltype>::fill(const soltype &val)
{
  std::fill(array.data(), array.data() + array.num_elements(), val);
}

template <class soltype>
void Integrator::History<soltype>::initialize_past(const soltype &val)
{
  for(int n = 0; n < static_cast<int>(array.shape()[PARTICLES]); ++n) {
    for(int t = array.index_bases()[TIMES]; t <= 0; ++t) {
      array[n][t][DERIV_0] = val;
    }
  }
}

#endif
