#include <boost/test/unit_test.hpp>
#include <iomanip>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/AIM/nearfield_interaction.h"
#include "interactions/direct_interaction.h"

BOOST_AUTO_TEST_SUITE(DIRECT_AIM_COMPARISON)

struct PARAMETERS {
  int interpolation_order, expansion_order, num_steps, num_dots;
  double c, dt, total_time, omega;

  Eigen::Array3d spacing;

  std::shared_ptr<DotVector> dots;
  std::shared_ptr<Integrator::History<Eigen::Vector2cd>> history;

  PARAMETERS()
      : interpolation_order(4),
        expansion_order(4),
        num_steps(1024),
        num_dots(2),

        c(1),
        dt(1),
        total_time(dt * num_steps),
        omega(M_PI / 10),

        spacing(Eigen::Array3d(1, 1, 1) * c * dt),

        dots(std::make_shared<DotVector>(num_dots)),
        history(std::make_shared<Integrator::History<Eigen::Vector2cd>>(
            num_dots, 10, num_steps))
  {
    history->fill(Eigen::Vector2cd::Zero());
    for(int d = 0; d < num_dots; ++d) {
      for(int i = -10; i < num_steps; ++i) {
        // Everybody radiates a Gaussian
        history->array[d][i][0](RHO_01) = src(i * dt);
      }
    }
  };

  double src(const double t) const
  {
    double arg = (t - total_time / 2.0) / (total_time / 6.0);
    return gaussian(arg);
  }
};

BOOST_FIXTURE_TEST_CASE(L2_ERROR, PARAMETERS)
{
  const Eigen::Vector3d zhat(0, 0, 1);

  dots->at(0) = QuantumDot(Eigen::Vector3d(0.2, 0.2, 0.2), zhat);
  dots->at(1) = QuantumDot(Eigen::Vector3d(0.7, 0.7, 0.7), zhat);

  auto greens_function = Propagation::RotatingFramePropagator(1, c, omega);

  DirectInteraction direct(dots, history, greens_function, interpolation_order,
                           c, dt);
  AIM::NearfieldInteraction nf(dots, history, greens_function,
                               interpolation_order, c, dt,
                               AIM::Grid(spacing, dots, expansion_order));

  for(int t = 0; t < num_steps; ++t) {
    std::cout << direct.evaluate(t).transpose() << " "
              << nf.evaluate(t).transpose() << std::endl;
  }
}
BOOST_AUTO_TEST_SUITE_END()  // DIRECT_AIM_COMPARISON
