#include "interaction_table.h"

using Vec3d = Eigen::Vector3d;

std::pair<int, double> compute_delay(const double);
Eigen::Matrix3d rhat_dyadic(const Vec3d &);
Eigen::Matrix3d nearfield_dyadic(const Vec3d &);
Eigen::Matrix3d midfield_dyadic(const Vec3d &);
Eigen::Matrix3d farfield_dyadic(const Vec3d &);

InteractionTable::InteractionTable(
    const int interp_order,
    const std::shared_ptr<const std::vector<QuantumDot>> &dots)
    : incident_interaction(dots->size()),
      history_interaction(dots->size()),
      interp_order(interp_order),
      num_interactions(dots->size() * (dots->size() - 1) / 2),
      dots(dots),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1])
{
  UniformLagrangeSet lagrange(interp_order, config.dt);
  for(size_t pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);

    Vec3d dr(separation((*dots)[src], (*dots)[obs]));

    std::pair<int, double> delay(
        compute_delay(dr.norm() / (config.c0 * config.dt)));

    floor_delays.at(pair_idx) = delay.first;
    lagrange.calculate_weights(delay.second);

    for(int i = 0; i <= interp_order; ++i) {
      coefficients[pair_idx][i] =
          dyadic_product((*dots)[obs], nearfield_dyadic(dr), (*dots)[src]) *
              lagrange.weights[0][i] +
          dyadic_product((*dots)[obs], midfield_dyadic(dr), (*dots)[src]) *
              lagrange.weights[1][i] +
          dyadic_product((*dots)[obs], farfield_dyadic(dr), (*dots)[src]) *
              lagrange.weights[2][i];
    }
  }
}

void InteractionTable::compute_interactions(
    const Pulse &pulse, const PredictorCorrector::HistoryArray &history,
    const int time_idx)
{
  compute_incident_interaction(pulse, time_idx * config.dt);
  compute_history_interaction(history, time_idx);
}

void InteractionTable::compute_incident_interaction(const Pulse &pulse,
                                                    const double time)
{
  for(size_t i = 0; i < dots->size(); ++i) {
    incident_interaction[i] =
        pulse((*dots)[i].position(), time).dot((*dots)[i].dipole());
  }
}

void InteractionTable::compute_history_interaction(
    const PredictorCorrector::HistoryArray &history, const int time_idx)
{
  for(size_t pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < history.index_bases()[1]) continue;
      history_interaction[src] +=
          polarization(history[obs][s - i][0]) * coefficients[pair_idx][i];
      history_interaction[obs] +=
          polarization(history[src][s - i][0]) * coefficients[pair_idx][i];
    }
  }
}

int InteractionTable::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<int, int> InteractionTable::idx2coord(const int idx)
{
  const int row = std::floor((std::sqrt(1 + 8 * idx) + 1) / 2.0);
  const int col = idx - row * (row - 1) / 2;

  return std::pair<int, int>(row, col);
}

std::pair<int, double> compute_delay(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

Eigen::Matrix3d rhat_dyadic(const Vec3d &dr)
{
  return dr * dr.transpose() / dr.squaredNorm();
}

Eigen::Matrix3d nearfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr));
  return dyad * std::pow(config.c0, 2) / std::pow(dr.norm(), 3);
}

Eigen::Matrix3d midfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr));
  return dyad * config.c0 / dr.squaredNorm();
}

Eigen::Matrix3d farfield_dyadic(const Vec3d &dr)
{
  Eigen::Matrix3d dyad = (Eigen::Matrix3d::Identity() - rhat_dyadic(dr));
  return dyad / dr.norm();
}
