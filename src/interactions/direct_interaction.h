#ifndef DIRECT_INTERACTION_H
#define DIRECT_INTERACTION_H

#include "history_interaction.h"

class DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::RotatingEFIE,
      const int,
      const double,
      const double);

  const ResultArray &evaluate(const int) final;

 private:
  Propagation::RotatingEFIE propagator;
  int num_interactions;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;

  void build_coefficient_table();

  static int coord2idx(int, int);
  static std::pair<int, int> idx2coord(const int);
};

#endif
