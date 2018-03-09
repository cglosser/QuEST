#ifndef NEARFIELD_H
#define NEARFIELD_H

#include "aim_base.h"

namespace AIM {
  class Nearfield;
}

class AIM::Nearfield final : public AimBase {
 public:
  Nearfield(std::shared_ptr<const DotVector>,
            std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
            const int,
            const int,
            const double,
            const double,
            const Grid &,
            std::shared_ptr<const Expansions::ExpansionTable>,
            Normalization::SpatialNorm,
            const boost::multi_array<double, 4> &,
            Projector::Projector_fn<cmplx>);

  ~Nearfield() = default;

 private:
  enum FIELD_AXIS_LABEL { STEPS, PAIRS, POINTSETS, EXPANSIONS, DIMS };
  std::array<int, 5> field_table_dims;
  std::vector<const_DotRange> mapping;
  std::vector<Grid::ipair_t> neighbors;

  spacetime::vector<cmplx> make_propagation_table() const override;
  void fill_source_table(const int) override;
  void propagate(const int) override;
  void fill_chebyshev_table(const int) override;
};

#endif
