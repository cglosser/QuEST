#include "history_interaction.h"

using Vec3d = Eigen::Vector3d;

HistoryInteraction::HistoryInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const Integrator::History<soltype>> &history,
    const std::shared_ptr<Propagation::FixedFramePropagator> &dyadic,
    const int interp_order,
    const double dt,
    const double c0)
    : Interaction(dots),
      history(history),
      dyadic(dyadic),
      interp_order(interp_order),
      num_interactions(dots->size() * (dots->size() - 1) / 2),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1]),
      dt(dt),
      c0(c0)
{
  build_coefficient_table();
}

void HistoryInteraction::build_coefficient_table()
{
  Interpolation::UniformLagrangeSet lagrange(interp_order);

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Vec3d dr(separation((*dots)[src], (*dots)[obs]));

    std::pair<int, double> delay(split_double(dr.norm() / (c0 * dt)));

    floor_delays[pair_idx] = delay.first;
    lagrange.calculate_weights(delay.second, dt);

    std::vector<double> interp_dyads(dyadic->coefficients(dr, lagrange));

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] = interp_dyads[i];
    }
  }
}

const Interaction::ResultArray &HistoryInteraction::evaluate(const int time_idx)
{
  for(unsigned int i = 0; i < results.size(); ++i)
    results[i] = Eigen::Vector3d(0, 0, 0);

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);
    const int s = time_idx - floor_delays[pair_idx];

    Vec3d dr(separation((*dots)[src], (*dots)[obs]));

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < history->array.index_bases()[1]) continue;

      results[src] +=
          -dr.cross(coefficients[pair_idx][i] * history->array[obs][s - i][0]);
      results[obs] +=
          -dr.cross(coefficients[pair_idx][i] * history->array[src][s - i][0]);
    }
  }

  return results;
}

int HistoryInteraction::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<int, int> HistoryInteraction::idx2coord(const int idx)
{
  const int row = std::floor((std::sqrt(1 + 8 * idx) + 1) / 2.0);
  const int col = idx - row * (row - 1) / 2;

  return std::pair<int, int>(row, col);
}

std::pair<int, double> HistoryInteraction::split_double(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}
