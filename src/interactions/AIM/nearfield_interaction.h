#pragma once

#include "interactions/AIM/grid.h"
#include "interactions/history_interaction.h"

namespace AIM {
  class NearfieldInteraction;
}

class AIM::NearfieldInteraction final : public HistoryInteraction {
 public:
  NearfieldInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::RotatingFramePropagator,
      const int,
      const double,
      const double,
      Grid);

  std::vector<std::pair<int, int>> build_pair_list() const;
  const ResultArray &evaluate(const int) final;

 private:
  Propagation::RotatingFramePropagator propagator;
  int num_interactions;
  Grid grid;
  std::vector<std::pair<int, int>> interaction_pairs;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;

  // boost multiarrays don't implement move constructors,
  // otherwise this would totally be up with build_pair_list.
  void build_coefficient_table();
};
