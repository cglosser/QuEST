#include "history_interaction.h"

using Vec3d = Eigen::Vector3d;

namespace History {
  std::pair<int, double> compute_delay(const double);
  Eigen::Matrix3d rhat_dyadic(const Vec3d &);
  Eigen::Matrix3d nearfield_dyadic(const Vec3d &);
  Eigen::Matrix3d midfield_dyadic(const Vec3d &);
  Eigen::Matrix3d farfield_dyadic(const Vec3d &);
}

HistoryInteraction::HistoryInteraction(
    const std::shared_ptr<const DotVector> &dots,
    const std::shared_ptr<const History::HistoryArray> &history,
    const int interp_order)
    : Interaction(dots),
      history(history),
      interp_order(interp_order),
      num_interactions(dots->size() * (dots->size() - 1) / 2),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1])
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

    std::pair<int, double> delay(
        History::compute_delay(dr.norm() / (config.c0 * config.dt)));

    floor_delays.at(pair_idx) = delay.first;
    lagrange.calculate_weights(delay.second);

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] =
          -config.mu0 / (4 * M_PI) *
          (dyadic_product((*dots)[obs], History::nearfield_dyadic(dr), (*dots)[src]) *
               lagrange.weights[0][i] +
           dyadic_product((*dots)[obs], History::midfield_dyadic(dr), (*dots)[src]) *
               lagrange.weights[1][i] / config.dt +
           dyadic_product((*dots)[obs], History::farfield_dyadic(dr), (*dots)[src]) *
               lagrange.weights[2][i] / std::pow(config.dt, 2));
    }
  }
}

void HistoryInteraction::evaluate(const int time_idx)
{
  std::fill(results.begin(), results.end(), 0);

  for(int pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < (*history).index_bases()[1]) continue;
      results[src] +=
          polarization((*history)[obs][s - i][0]) * coefficients[pair_idx][i];
      results[obs] +=
          polarization((*history)[src][s - i][0]) * coefficients[pair_idx][i];
    }
  }
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

std::pair<int, double> History::compute_delay(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

Eigen::Matrix3d History::rhat_dyadic(const Vec3d &dr)
{
  return dr * dr.transpose() / dr.squaredNorm();
}

Eigen::Matrix3d History::nearfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - 3 * History::rhat_dyadic(dr));
  return dyad * std::pow(config.c0, 2) / std::pow(dr.norm(), 3);
}

Eigen::Matrix3d History::midfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - 3 * History::rhat_dyadic(dr));
  return dyad * config.c0 / dr.squaredNorm();
}

Eigen::Matrix3d History::farfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - History::rhat_dyadic(dr));
  return dyad / dr.norm();
}
