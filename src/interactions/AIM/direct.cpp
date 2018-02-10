#include "interactions/AIM/direct.h"

AIM::DirectInteraction::DirectInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::Kernel<cmplx> &kernel,
    const int interp_order,
    const double c0,
    const double dt,
    const Grid &grid)
    : HistoryInteraction(dots, history, interp_order, c0, dt),
      interaction_pairs(make_pair_list(grid)),
      shape({static_cast<int>(interaction_pairs.size()), interp_order + 1}),
      floor_delays(shape[0]),
      coefficients(shape)
{
  build_coefficient_table(kernel);
}

const Interaction::ResultArray &AIM::DirectInteraction::evaluate(
    const int time_idx)
{
  results.setZero();

  for(int pair_idx = 0; pair_idx < shape[0]; ++pair_idx) {
    const auto &pair = interaction_pairs[pair_idx];
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i < shape[1]; ++i) {
      if(s - i < history->array_.index_bases()[1]) continue;

      constexpr int RHO_01 = 1;

      results[pair.first] += (history->array_[pair.second][s - i][0])[RHO_01] *
                             coefficients[pair_idx][i];
      results[pair.second] += (history->array_[pair.first][s - i][0])[RHO_01] *
                              coefficients[pair_idx][i];
    }
  }

  return results;
}

std::vector<std::pair<int, int>> AIM::DirectInteraction::make_pair_list(
    const Grid &grid) const
{
  std::vector<std::pair<int, int>> particle_pairs;
  std::vector<DotRange> mapping{grid.box_contents_map()};

  for(const auto &p : grid.nearfield_pairs(1)) {
    DotVector::const_iterator begin1, end1, begin2, end2;

    std::tie(begin1, end1) = mapping[p.first];
    std::tie(begin2, end2) = mapping[p.second];

    for(auto dot1 = begin1; dot1 != end1; ++dot1) {
      auto idx1{std::distance(dots->begin(), dot1)};
      for(auto dot2 = begin2; dot2 != end2; ++dot2) {
        auto idx2{std::distance(dots->begin(), dot2)};

        particle_pairs.emplace_back(idx1, idx2);
      }
    }
  }

  return particle_pairs;
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
