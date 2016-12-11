#include "interpolation_table.h"

using Vec3d = Eigen::Vector3d;

std::pair<int, double> compute_delay(const double);
Eigen::Matrix3d rhat_dyadic(const Vec3d &);
Eigen::Matrix3d nearfield_dyadic(const Vec3d &);
Eigen::Matrix3d midfield_dyadic(const Vec3d &);
Eigen::Matrix3d farfield_dyadic(const Vec3d &);

InterpolationTable::InterpolationTable(
    const int n, std::shared_ptr<const std::vector<QuantumDot>> qdots)
    : convolution(qdots->size()),
      interp_order(n),
      num_interactions(qdots->size() * (qdots->size() - 1) / 2),
      dots(qdots),
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

void InterpolationTable::compute_interactions(const Pulse &pulse,
                                              const HistoryArray &history,
                                              const int time_idx)
{
  compute_incident_interaction(pulse, time_idx * config.dt);
  convolve_currents(history, time_idx);
}

void InterpolationTable::compute_incident_interaction(const Pulse &pulse,
                                                      const double time)
{
  for(size_t i = 0; i < dots->size(); ++i) {
    convolution[i] =
        pulse((*dots)[i].position(), time).dot((*dots)[i].dipole());
  }
}

void InterpolationTable::convolve_currents(const HistoryArray &history,
                                           const int time_idx)
{
  for(size_t pair_idx = 0; pair_idx < num_interactions; ++pair_idx) {
    int src, obs;
    std::tie(src, obs) = idx2coord(pair_idx);
    const int s = time_idx - floor_delays[pair_idx];

    for(int i = 0; i <= interp_order; ++i) {
      if(s - i < history.index_bases()[1]) continue;
      convolution[src] +=
          polarization(history[obs][s - i][0]) * coefficients[pair_idx][i];
      convolution[obs] +=
          polarization(history[src][s - i][0]) * coefficients[pair_idx][i];
    }
  }
}

int InterpolationTable::coord2idx(int row, int col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<int, int> InterpolationTable::idx2coord(const int idx)
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
