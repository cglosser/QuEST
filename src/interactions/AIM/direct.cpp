#include "interactions/AIM/direct.h"

AIM::DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const int border,
    const double c0,
    const double dt,
    const Grid &grid)
    : HistoryInteraction(dots, history, interp_order, c0, dt),
      interaction_pairs(grid.nearfield_point_pairs(border, *dots)),
      shape({{static_cast<int>(interaction_pairs.size()), interp_order + 1}}),
      floor_delays(shape[0]),
      coefficients(shape)
{
  build_coefficient_table(kernel);
}

const InteractionBase::ResultArray &AIM::DirectInteraction::evaluate(
    const int time_idx)
{
  results.setZero();

  for(int pair_idx = 0; pair_idx < shape[0]; ++pair_idx) {
    const auto &pair = interaction_pairs[pair_idx];

    for(int i = 0; i < shape[1]; ++i) {
      const int s =
          std::max(time_idx - floor_delays[pair_idx] - i,
                   static_cast<int>(history->array_.index_bases()[1]));

      constexpr int RHO_01 = 1;

      results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                             coefficients[pair_idx][i];
      results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                              coefficients[pair_idx][i];
    }
  }

  return results;
}

void AIM::DirectInteraction::build_coefficient_table(
    Propagation::Kernel<cmplx> &kernel)
{
  Interpolation::UniformLagrangeSet lagrange(shape[1] - 1);

  for(int pair_idx = 0; pair_idx < shape[0]; ++pair_idx) {
    const std::pair<const int, const int> &pair = interaction_pairs[pair_idx];

    Eigen::Vector3d dr(separation((*dots)[pair.first], (*dots)[pair.second]));
    auto delay = split_double(dr.norm() / (c0 * dt));

    floor_delays[pair_idx] = delay.first;
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        kernel.coefficients(dr, lagrange));

    for(int i = 0; i < shape[1]; ++i) {
      coefficients[pair_idx][i] = (*dots)[pair.second].dipole().dot(
          interp_dyads[i] * (*dots)[pair.first].dipole());
    }
  }
}
