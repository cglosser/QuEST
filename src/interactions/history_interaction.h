#ifndef PULSE_INTERACTION_H
#define PULSE_INTERACTION_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include "../history.h"
#include "../configuration.h"
#include "../lagrange_set.h"
#include "../quantum_dot.h"
#include "interaction.h"

class HistoryInteraction : public Interaction {
 public:
  HistoryInteraction(const std::shared_ptr<const DotVector> &,
                     const std::shared_ptr<const History::HistoryArray> &,
                     const int);
  void evaluate(const int);

 private:
  std::shared_ptr<const History::HistoryArray> history;
  int interp_order, num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<double, 2> coefficients;

  void build_coefficient_table();

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
