#include "nearfield.h"

AIM::Nearfield::Nearfield(
    const std::shared_ptr<const DotVector> dots,
    const std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    const int interp_order,
    const int border,
    const double c0,
    const double dt,
    std::shared_ptr<const Grid> grid,
    std::shared_ptr<const Expansions::ExpansionTable> expansion_table,
    Expansions::ExpansionFunction expansion_function,
    Normalization::SpatialNorm normalization)
    : AimBase(dots,
              history,
              interp_order,
              c0,
              dt,
              grid,
              expansion_table,
              expansion_function,
              normalization),
      interaction_pairs_(make_pair_list(border)),
      shape_({{static_cast<int>(interaction_pairs_.size()),
               grid->max_transit_steps(c0, dt) + interp_order}}),
      coefficients_{build_coefficient_table(interp_order)}
{
}

const InteractionBase::ResultArray &AIM::Nearfield::evaluate(const int time_idx)
{
  results.setZero();

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = interaction_pairs_[pair_idx];

    for(int i = 0; i < shape_[1]; ++i) {
      const int s = std::max(
          time_idx - i, static_cast<int>(history->array_.index_bases()[1]));

      constexpr int RHO_01 = 1;

      results[pair.first] += (history->array_[pair.second][s][0])[RHO_01] *
                             coefficients_[pair_idx][i];
      results[pair.second] += (history->array_[pair.first][s][0])[RHO_01] *
                              coefficients_[pair_idx][i];
    }
  }

  return results;
}

std::vector<std::pair<int, int>> AIM::Nearfield::make_pair_list(
    const int border) const
{
  std::vector<std::pair<int, int>> particle_pairs;
  auto mapping = grid->box_contents_map(*dots);

  for(const auto &p : grid->nearfield_pairs(border, *dots)) {
    DotVector::const_iterator begin1, end1, begin2, end2;

    std::tie(begin1, end1) = mapping[p.first];
    std::tie(begin2, end2) = mapping[p.second];

    for(auto dot1 = begin1; dot1 != end1; ++dot1) {
      auto idx1{std::distance(dots->begin(), dot1)};
      for(auto dot2 = begin2; dot2 != end2; ++dot2) {
        auto idx2{std::distance(dots->begin(), dot2)};

        if(idx1 <= idx2) particle_pairs.emplace_back(idx1, idx2);
      }
    }
  }

  return particle_pairs;
}

boost::multi_array<cmplx, 2> AIM::Nearfield::build_coefficient_table(
    const int interp_order) const
{
  boost::multi_array<cmplx, 2> coefficients(shape_);
  std::fill(coefficients.data(),
            coefficients.data() + coefficients.num_elements(), cmplx(0, 0));

  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int pair_idx = 0; pair_idx < shape_[0]; ++pair_idx) {
    const auto &pair = interaction_pairs_[pair_idx];

    const double innerprod =
        ((pair.first == pair.second) ? 0.5 : 1) *
        (*dots)[pair.second].dipole().dot((*dots)[pair.first].dipole());

    for(size_t e1 = 0; e1 < expansion_table->shape()[1]; ++e1) {
      for(size_t e2 = 0; e2 < expansion_table->shape()[1]; ++e2) {
        int idx1 = (*expansion_table)[pair.first][e1].index,
            idx2 = (*expansion_table)[pair.second][e2].index;

        if(idx1 == idx2) continue;

        Eigen::Vector3d dr(grid->spatial_coord_of_box(idx2) -
                           grid->spatial_coord_of_box(idx1));

        const double arg = dr.norm() / (c0 * dt);
        const auto split_arg = split_double(arg);

        for(int t = 1; t < shape[1]; ++t) {
          const auto polynomial_idx = static_cast<int>(ceil(t - arg));
          if(0 <= polynomial_idx && polynomial_idx <= interp_order) {
            lagrange.evaluate_derivative_table_at_x(split_arg.second, dt);
            cmplx matrix_element =
                lagrange.evaluations[0][polynomial_idx] * normalization(dr);

            coefficients[pair_idx][t] +=
                innerprod * (*expansion_table)[pair.second][e2].d0 *
                matrix_element * (*expansion_table)[pair.first][e1].d0;
          }
        }
      }
    }
  }

  return coefficients;
}
