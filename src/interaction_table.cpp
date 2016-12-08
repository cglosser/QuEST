#include "interaction_table.h"

using Vec3d = Eigen::Vector3d;

std::pair<int, double> compute_delay(const double);
Eigen::Matrix3d rhat_dyadic(const Vec3d &);
Eigen::Matrix3d nearfield_dyadic(const Vec3d &);
Eigen::Matrix3d midfield_dyadic(const Vec3d &);
Eigen::Matrix3d farfield_dyadic(const Vec3d &);

InteractionTable::InteractionTable(const int n, DotTable qdots)
    : convolution(qdots.size()),
      interp_order(n),
      num_interactions(qdots.size() * (qdots.size() - 1) / 2),
      dots(qdots),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1])
{
  UniformLagrangeSet lagrange(interp_order, config.dt);
  for(size_t src = 0; src < qdots.size() - 1; ++src) {
    for(size_t obs = src + 1; obs < qdots.size(); ++obs) {
      int idx = coord2idx(src, obs);

      Vec3d dr(separation(qdots[src], qdots[obs]));

      std::pair<int, double> delay(
          compute_delay(dr.norm() / (config.c0 * config.dt)));

      floor_delays.at(idx) = delay.first;
      lagrange.calculate_weights(delay.second);

      for(int i = 0; i <= interp_order; ++i) {
        coefficients[idx][i] =
            dyadic_product(qdots[obs], nearfield_dyadic(dr), qdots[src]) *
                lagrange.weights[0][i] +
            dyadic_product(qdots[obs], midfield_dyadic(dr), qdots[src]) *
                lagrange.weights[1][i] +
            dyadic_product(qdots[obs], farfield_dyadic(dr), qdots[src]) *
                lagrange.weights[2][i];
      }
    }
  }
}

void InteractionTable::convolve_currents(const HistoryArray &history,
                                         const int time_idx)
{
  std::fill(convolution.begin(), convolution.end(), 0);

  for(size_t src = 0; src < dots.size() - 1; ++src) {
    for(size_t obs = src + 1; obs < dots.size(); ++obs) {
      const int idx = coord2idx(src, obs);
      const int s = time_idx - floor_delays[idx];

      for(int i = 0; i <= interp_order; ++i) {
        if(s - i < history.index_bases()[1]) continue;
        convolution[src] +=
            polarization(history[obs][s - i][0]) * coefficients[idx][i];
        convolution[obs] +=
            polarization(history[src][s - i][0]) * coefficients[idx][i];
      }
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
