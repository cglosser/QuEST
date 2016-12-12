#ifndef INTERACTION_TABLE_H
#define INTERACTION_TABLE_H

#include <Eigen/Dense>
#include <algorithm>
#include <boost/multi_array.hpp>
#include <cmath>
#include <tuple>
#include <vector>

#include "configuration.h"
#include "integrator.h"
#include "lagrange_set.h"
#include "pulse.h"
#include "quantum_dot.h"

class InteractionTable {
 public:
  InteractionTable(const int, const std::shared_ptr<const DotVector> &,
                   const std::shared_ptr<const Pulse> &);
  std::vector<double> incident_interaction, history_interaction;

  void compute_interactions(const PredictorCorrector::HistoryArray &,
                            const int);
  void compute_history_interaction(const PredictorCorrector::HistoryArray &,
                                   const int);

  double result(const int i)
  {
    return incident_interaction[i] + history_interaction[i];
  }

 private:
  int interp_order, num_interactions;
  std::shared_ptr<const DotVector> dots;
  std::shared_ptr<const Pulse> pulse;
  std::vector<int> floor_delays;
  boost::multi_array<double, 2> coefficients;

  void compute_incident_interaction(const double);
  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
