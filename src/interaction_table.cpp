#include "interaction_table.h"

std::pair<int, double> compute_delay(const double);
Eigen::Matrix3d rhat_dyadic(const Eigen::Vector3d &);
double nearfield_dyadic(const Eigen::Vector3d &, const Eigen::Vector3d &,
                        const Eigen::Vector3d &);
double midfield_dyadic(const Eigen::Vector3d &, const Eigen::Vector3d &,
                       const Eigen::Vector3d &);
double farfield_dyadic(const Eigen::Vector3d &, const Eigen::Vector3d &,
                       const Eigen::Vector3d &);

InteractionTable::InteractionTable(const int n, const std::vector<QuantumDot> &dots)
    : interp_order(n),
      num_dots(dots.size()),
      num_interactions(dots.size()*(dots.size() - 1)/2),
      coefficients(boost::extents[num_interactions][interp_order + 1])
{
  UniformLagrangeSet lagrange(interp_order);
  for(size_t src = 0; src < dots.size() - 1; ++src) {
    for(size_t obs = src + 1; obs < dots.size(); ++obs) {
      Eigen::Vector3d dr(dots[obs].pos - dots[src].pos);
      size_t idx = coord2idx(src, obs);

      std::pair<int, double>
        delay(compute_delay(dr.norm()/(config.c0*config.dt)));

      lagrange.calculate_weights(delay.second);

      for(int i = 0; i <= interp_order; ++i) {
        coefficients[idx][i] =
            nearfield_dyadic(dr, dots[src].dipole, dots[obs].dipole)*lagrange.weights[0][i] +
            midfield_dyadic(dr, dots[src].dipole, dots[obs].dipole)*lagrange.weights[1][i] +
            farfield_dyadic(dr, dots[src].dipole, dots[obs].dipole)*lagrange.weights[2][i];
      }
    }
  }
}

size_t InteractionTable::coord2idx(size_t row, size_t col)
{
  assert(row != col);
  if(col < row) std::swap(row, col);

  return num_interactions -
         (num_dots - row) * (num_dots - row - 1) / 2 + col - row - 1;
}

std::pair<size_t, size_t> InteractionTable::idx2coord(const size_t idx)
{
  const size_t r =
      num_dots - 2 -
      std::floor(std::sqrt(-8 * idx + 8 * num_interactions - 7) / 2 -
                 0.5);
  const size_t c = idx + r + 1 - num_dots * (num_dots - 1) / 2 +
                   (num_dots - r) * (num_dots - r - 1) / 2;

  return std::pair<size_t, size_t>(r, c);
}

std::pair<int, double> compute_delay(const double delay)
{
  std::pair<int, double> result;

  double idelay;
  result.second = std::modf(delay, &idelay);
  result.first = static_cast<int>(idelay);

  return result;
}

Eigen::Matrix3d rhat_dyadic(const Eigen::Vector3d &dr)
{
  return dr * dr.transpose() / dr.squaredNorm();
}

double nearfield_dyadic(const Eigen::Vector3d &dr, const Eigen::Vector3d &src,
                        const Eigen::Vector3d &obs)
{
  const double dyad = obs.transpose() *
                      (Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr)) * src;
  return dyad * std::pow(config.c0, 2) / std::pow(dr.norm(), 3);
}

double midfield_dyadic(const Eigen::Vector3d &dr, const Eigen::Vector3d &src,
                       const Eigen::Vector3d &obs)
{
  const double dyad = obs.transpose() *
                      (Eigen::Matrix3d::Identity() - 3 * rhat_dyadic(dr)) * src;
  return dyad * config.c0 / dr.squaredNorm();
}

double farfield_dyadic(const Eigen::Vector3d &dr, const Eigen::Vector3d &src,
                       const Eigen::Vector3d &obs)
{
  const double dyad =
      obs.transpose() * (Eigen::Matrix3d::Identity() - rhat_dyadic(dr)) * src;
  return dyad / dr.norm();
}
