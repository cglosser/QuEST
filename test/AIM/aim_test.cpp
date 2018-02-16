#include <boost/test/unit_test.hpp>
#include <cmath>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

template <int num_dots>
struct PARAMETERS {
  int interpolation_order, expansion_order, border, num_steps;
  double c, dt, total_time;

  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(1),
        border(1),
        num_steps(1024),

        c(1),
        dt(1),
        total_time(1024),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>(num_dots)),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps))
  {
    std::cout.precision(17);

    history->fill(Eigen::Vector2cd::Zero());
    for(int t = -10; t < num_steps; ++t) {
      history->array_[0][t][0](RHO_01) = src(t * dt);
    }
  };

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return Math::gaussian(arg);
  }
};

BOOST_AUTO_TEST_SUITE(IDENTITY_KERNEL)

BOOST_AUTO_TEST_SUITE_END()  // IDENTITY_KERNEL

BOOST_AUTO_TEST_SUITE_END()  // AIM
