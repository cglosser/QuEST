#pragma once

#include "interactions/AIM/grid.h"
#include "interactions/history_interaction.h"

namespace AIM {
  class DirectInteraction;
}

class AIM::DirectInteraction final : public HistoryInteraction {
 public:
  DirectInteraction(
      std::shared_ptr<const DotVector>,
      std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
      Propagation::Kernel<cmplx> &,
      const int,
      const int,
      const double,
      const double,
      const AIM::Grid &);

  const ResultArray &evaluate(const int) final;

 private:
  std::vector<std::pair<int, int>> interaction_pairs;
  std::array<int, 2> shape;
  std::vector<int> floor_delays;
  boost::multi_array<cmplx, 2> coefficients;

  void build_coefficient_table(Propagation::Kernel<cmplx> &);
};
