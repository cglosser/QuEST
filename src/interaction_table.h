#ifndef INTERACTION_TABLE_H
#define INTERACTION_TABLE_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <cmath>
#include <vector>

#include "configuration.h"
#include "lagrange_set.h"
#include "quantum_dot.h"

class InteractionTable {
 public:
  typedef boost::multi_array<matrix_elements, 3,
                             Eigen::aligned_allocator<matrix_elements>>
      HistoryArray;

  InteractionTable(const int, DotTable);

  std::vector<double> convolution;
  void convolve_currents(const HistoryArray &, const int);

 private:
  int interp_order;
  int num_interactions;
  DotTable dots;

  std::vector<int> floor_delays;
  boost::multi_array<double, 2> coefficients;

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
