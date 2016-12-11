#ifndef INTERPOLATION_TABLE_H
#define INTERPOLATION_TABLE_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <cmath>
#include <tuple>
#include <vector>

#include "configuration.h"
#include "lagrange_set.h"
#include "pulse.h"
#include "quantum_dot.h"

class InterpolationTable {
 public:
  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  InterpolationTable(const int, std::vector<QuantumDot>);
  std::vector<double> convolution;

  void compute_interactions(const Pulse &, const HistoryArray &, const int);

 private:
  int interp_order;
  int num_interactions;
  std::vector<QuantumDot> dots;

  std::vector<int> floor_delays;
  boost::multi_array<double, 2> coefficients;

  void compute_incident_interaction(const Pulse &, const double);
  void convolve_currents(const HistoryArray &, const int);

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
