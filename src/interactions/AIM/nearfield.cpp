#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const double c0,
    const double dt,
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const Expansions::ExpansionTable> expansion_table,
    Expansions::ExpansionFunction expansion_function,
    Normalization::SpatialNorm normalization,
    std::shared_ptr<const std::vector<Grid::ipair_t>> interaction_pairs,
    const double omega)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
              grid,
              expansion_table,
              expansion_function,
              normalization),
      omega_{omega},
      interaction_pairs_{std::move(interaction_pairs)},
      shape_({{static_cast<int>(interaction_pairs_->size()),
               grid->max_transit_steps(c0, dt) + 6, 2}}),
      coefficients_{coefficient_table()}
{
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate(const int time_idx)
{
  results.setZero();
  constexpr int RHO_01 = 1;

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];

    for(int t = 0; t < shape_[1]; ++t) {
      const int s = std::max(
          time_idx - t, static_cast<int>(history->array_.index_bases()[1]));

      results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                             coefficients_[pair_idx][t][0];

      results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                              coefficients_[pair_idx][t][1];
    }
  }

  return results;
}

boost::multi_array<cmplx, 3> AIM::Nearfield::coefficient_table() const
{
  boost::multi_array<cmplx, 3> coefficients(shape_);
  std::fill(coefficients.data(),
            coefficients.data() + coefficients.num_elements(), cmplx(0, 0));

  Interpolation::DerivFive lagrange(dt);
  const double c_sq = std::pow(c0, 2);

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = (*interaction_pairs_)[pair_idx];
    const auto &dot0 = (*dots)[pair.first], &dot1 = (*dots)[pair.second];

    for(size_t i = 0; i < expansion_table->shape()[1]; ++i) {
      const auto &e0 = (*expansion_table)[pair.first][i];

      for(size_t j = 0; j < expansion_table->shape()[1]; ++j) {
        const auto &e1 = (*expansion_table)[pair.second][j];

        if(e0.index == e1.index) continue;

        Eigen::Vector3d dr(grid->spatial_coord_of_box(e1.index) -
                           grid->spatial_coord_of_box(e0.index));

        cmplx matrix_element = normalization(dr);

        const double arg = dr.norm() / (c0 * dt);
        const auto split_arg = split_double(arg);

        lagrange.evaluate_derivative_table_at_x(split_arg.second);

        for(int poly = 0; poly < lagrange.order() + 1; ++poly) {
          // const cmplx time =
          //(lagrange.evaluations[2][poly] +
          // 2.0 * iu * omega_ * lagrange.evaluations[1][poly] -
          // std::pow(omega_, 2) * lagrange.evaluations[0][poly]) *
          // innerprod;

          coefficients[pair_idx][split_arg.first + poly][0] +=
              matrix_element * dot0.dipole().dot(e1.d0 * dot1.dipole()) *
              lagrange.evaluations[0][poly];

          if(pair.first == pair.second) continue;

          coefficients[pair_idx][split_arg.first + poly][1] +=
              matrix_element * dot1.dipole().dot(e0.d0 * dot0.dipole()) *
              lagrange.evaluations[0][poly];

          // coefficients[pair_idx][split_arg.first + poly][0] +=
          //-matrix_element *
          //(time - c_sq * del_sq[0] * lagrange.evaluations[0][poly]);

          // if(pair.first == pair.second) continue;

          // coefficients[pair_idx][split_arg.first + poly][1] +=
          //-matrix_element *
          //(time - c_sq * del_sq[1] * lagrange.evaluations[0][poly]);
        }
      }
    }
  }

  return coefficients;
}
