#include "interaction_table.h"

using Vec3d = Eigen::Vector3d;

std::pair<int, double> compute_delay(const double);
Eigen::Matrix3d rhat_dyadic(const Vec3d &);
Eigen::Matrix3d nearfield_dyadic(const Vec3d &);
Eigen::Matrix3d midfield_dyadic(const Vec3d &);
Eigen::Matrix3d farfield_dyadic(const Vec3d &);

InteractionTable::InteractionTable(const int n,
                                   const std::vector<QuantumDot> &dots)
    : interp_order(n),
      num_interactions(dots.size() * (dots.size() - 1) / 2),
      floor_delays(num_interactions),
      coefficients(boost::extents[num_interactions][interp_order + 1])
{
  UniformLagrangeSet lagrange(interp_order);
  for(size_t src = 0; src < dots.size() - 1; ++src) {
    for(size_t obs = src + 1; obs < dots.size(); ++obs) {
      size_t idx = coord2idx(src, obs);

      Vec3d dr(separation(dots[src], dots[obs]));

      std::pair<int, double> delay(
          compute_delay(dr.norm() / (config.c0 * config.dt)));

      floor_delays.at(idx) = delay.first;
      lagrange.calculate_weights(delay.second);

      for(int i = 0; i <= interp_order; ++i) {
        coefficients[idx][i] =
            dyadic_product(dots[obs], nearfield_dyadic(dr), dots[src]) *
                lagrange.weights[0][i] +
            dyadic_product(dots[obs], midfield_dyadic(dr), dots[src]) *
                lagrange.weights[1][i]/config.dt +
            dyadic_product(dots[obs], farfield_dyadic(dr), dots[src]) *
                lagrange.weights[2][i]/std::pow(config.dt, 2);
      }
    }
  }
}

size_t InteractionTable::coord2idx(size_t row, size_t col)
{
  assert(row != col);
  if(col > row) std::swap(row, col);

  return row * (row - 1) / 2 + col;
}

std::pair<size_t, size_t> InteractionTable::idx2coord(const size_t idx)
{
  const size_t row = std::floor((std::sqrt(1 + 8 * idx) + 1) / 2.0);
  const size_t col = idx - row * (row - 1) / 2;

  return std::pair<size_t, size_t>(row, col);
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
