#include "nearfield_interaction.h"

AIM::NearfieldInteraction::NearfieldInteraction(
    std::shared_ptr<const DotVector> dots,
    std::shared_ptr<const Integrator::History<Eigen::Vector2cd>> history,
    Propagation::RotatingFramePropagator propagator,
    const int interp_order,
    const double c0,
    const double dt,
    Grid grid)
    : HistoryInteraction(
          std::move(dots), std::move(history), interp_order, c0, dt),
      propagator(std::move(propagator)),
      grid(std::move(grid)),
      interaction_pairs(build_pair_list()),
      floor_delays(interaction_pairs.size()),
      coefficients(boost::extents[interaction_pairs.size()][interp_order + 1])
{
  build_coefficient_table();
}

std::vector<std::pair<int, int>> AIM::NearfieldInteraction::build_pair_list()
    const
{
  std::vector<std::pair<int, int>> pairs;
  const auto box_contents = grid.box_contents_map();

  const auto get_dot_idx = [&](const DotVector::const_iterator &d) {
    return std::distance<DotVector::const_iterator>(dots->begin(), d);
  };

  for(auto box1 = box_contents.begin(); box1 != box_contents.end(); ++box1) {
    // start iter = end iter, thus empty box
    if(box1->first == box1->second) continue;

    for(auto src_dot = box1->first; src_dot != box1->second - 1; ++src_dot) {
      for(auto obs_dot = src_dot + 1; obs_dot != box1->second; ++obs_dot) {
        pairs.push_back({get_dot_idx(src_dot), get_dot_idx(obs_dot)});
      }
    }

    for(auto box2 = box1 + 1; box2 < box_contents.end(); ++box2) {
      if(box2->first == box2->second) continue;

      // box1 and box2 are *iterators over box_contents.* box_contents itself
      // contains *pairs of iterators* denoting the range of dots that exist
      // within that box/are associated with that gridpoint (assuming the
      // DotVector is sorted by that criterion). If the box is not empty (i.e.
      // *->first != *->second, checked above), then any of the particles within
      // the box can be used to determine the box's location and it's convenient
      // to choose the first particle)
      bool is_in_nearfield = grid.is_nearfield_pair(box1->first->position(),
                                                    box2->first->position());
      if(!is_in_nearfield) continue;

      for(auto src_dot = box1->first; src_dot != box1->second; ++src_dot) {
        for(auto obs_dot = box2->first; obs_dot != box2->second; ++obs_dot) {
          pairs.push_back({get_dot_idx(src_dot), get_dot_idx(obs_dot)});
        }
      }
    }
  }

  pairs.shrink_to_fit();
  return pairs;
}

void AIM::NearfieldInteraction::build_coefficient_table()
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(auto pair_idx = 0u; pair_idx < interaction_pairs.size(); ++pair_idx) {
    std::pair<int, int> &pair = interaction_pairs[pair_idx];

    Eigen::Vector3d dr(separation((*dots)[pair.first], (*dots)[pair.second]));
    auto delay = split_double(dr.norm() / (c0 * dt));

    floor_delays[pair_idx] = delay.first;
    lagrange.evaluate_derivative_table_at_x(delay.second, dt);

    std::vector<Eigen::Matrix3cd> interp_dyads(
        propagator.coefficients(dr, lagrange));

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] = (*dots)[pair.second].dipole().dot(
          interp_dyads[i] * (*dots)[pair.first].dipole());
    }
  }
}

const Interaction::ResultArray &AIM::NearfieldInteraction::evaluate(
    const int time_idx)
{
  results.setZero();

  for(auto pair_idx = 0u; pair_idx < interaction_pairs.size(); ++pair_idx) {
    const auto &pair = interaction_pairs[pair_idx];
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < history->array.index_bases()[1]) continue;

      constexpr int RHO_01 = 1;

      results[pair.first] += (history->array[pair.second][s - i][0])[RHO_01] *
                             coefficients[pair_idx][i];
      results[pair.second] += (history->array[pair.first][s - i][0])[RHO_01] *
                              coefficients[pair_idx][i];
    }
  }

  return results;
}
