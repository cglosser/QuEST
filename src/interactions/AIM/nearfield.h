#ifndef NEARFIELD_H
#define NEARFIELD_H

#include "aim_base.h"

namespace AIM {
  class Nearfield;
}

class AIM::Nearfield final : public AimBase {
 public:
  Nearfield(const std::shared_ptr<const DotVector>,
            const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
            const int,
            const double,
            const double,
            std::shared_ptr<const Grid>,
            std::shared_ptr<const Expansions::ExpansionTable>,
            Expansions::ExpansionFunction,
            Normalization::SpatialNorm,
            std::shared_ptr<const std::vector<Grid::ipair_t>>);
  ~Nearfield() = default;

  const ResultArray &evaluate(const int) final;

 private:
  std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs_;
  std::array<int, 2> shape_;
  boost::multi_array<cmplx, 2> coefficients_;

  boost::multi_array<cmplx, 2> coefficient_table(const int) const;
};

#endif
