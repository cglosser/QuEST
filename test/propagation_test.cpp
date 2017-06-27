#include <boost/test/unit_test.hpp>

#include "../src/propagation/history_evaluator.h"

BOOST_AUTO_TEST_SUITE(propagation)

struct System {
  double dt;
  int num_steps, interp_order;
  Integrator::History<double> history;
  std::shared_ptr<DotVector> dots;
  std::shared_ptr<GreenFunction::Dyadic> dyadic_gf;

  System()
      : dt(1e-2),
        num_steps(2000),
        interp_order(4),
        history(2, 20, num_steps),
        dots(std::make_shared<DotVector>(2)),
        dyadic_gf(std::make_shared<GreenFunction::Dyadic>(4 * M_PI, 1, 1))
  {
    history.fill(0);
    dots->at(0) =
        QuantumDot(Eigen::Vector3d(0, 0, 0), 2278.9013,
                   std::make_pair(10.0, 20.0), Eigen::Vector3d(0, 0, 1));
    dots->at(1) =
        QuantumDot(Eigen::Vector3d(1, 0, 0), 2278.9013,
                   std::make_pair(10.0, 20.0), Eigen::Vector3d(0, 0, 1));
  };
};

BOOST_FIXTURE_TEST_CASE(dyad, System)
{
  Eigen::Vector3d dr(dots->at(1).position() - dots->at(0).position());
  //UniformLagrangeSet uls(dr.norm()/dt
}

BOOST_AUTO_TEST_SUITE_END()
