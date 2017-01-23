#ifndef HISTORY_INTERACTION_H
#define HISTORY_INTERACTION_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

#include "../configuration.h"
#include "../history.h"
#include "../lagrange_set.h"
#include "../quantum_dot.h"
#include "green_function.h"
#include "interaction.h"

class HistoryInteraction : public Interaction {
 public:
  HistoryInteraction(const std::shared_ptr<const DotVector> &,
                     const std::shared_ptr<const History::HistoryArray> &,
                     const std::shared_ptr<GreenFunction::Dyadic> &,
                     const int);

  virtual const ResultArray &evaluate(const int);

 private:
  std::shared_ptr<const History::HistoryArray> history;
  std::shared_ptr<GreenFunction::Dyadic> dyadic;
  int interp_order, num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;

  void build_coefficient_table();

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
  static std::pair<int, double> split_double(const double);
};

#endif
