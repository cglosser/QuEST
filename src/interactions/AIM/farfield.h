#ifndef FARFIELD_H
#define FARFIELD_H

#include "aim_base.h"
#include "fourier.h"

namespace AIM {
  class Farfield;
}

class AIM::Farfield final : public AimBase {
 public:
  Farfield(const std::shared_ptr<const DotVector>,
           const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>>,
           const int,
           const double,
           const double,
           const Grid,
           const Expansions::ExpansionTable &,
           Expansions::ExpansionFunction,
           normalization::SpatialNorm);

  ~Farfield() = default;

 private:
  TransformPair spatial_vector_transforms;

  void fill_gmatrix_table(spacetime::vector<cmplx> &) const;
  TransformPair spatial_fft_plans();

  spacetime::vector<cmplx> make_propagation_table() override;
  void fill_source_table(const int) override;
  void propagate(const int) override;
  void fill_results_table(const int) override;
};

#endif
