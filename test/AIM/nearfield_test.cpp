#include <boost/test/unit_test.hpp>

#include "interactions/AIM/nearfield_interaction.h"

BOOST_AUTO_TEST_SUITE(NEARFIELD)

BOOST_AUTO_TEST_CASE(CHAIN)
{
  auto dots = std::make_shared<DotVector>();
  for(int i = 0; i < 20; ++i) {
    dots->push_back(QuantumDot(Eigen::Vector3d(0, 0, 0.1 + 0.7 * i)));
  }
  Propagation::RotatingFramePropagator g(1, 1, 0);

  AIM::Grid grid(Eigen::Vector3d(1, 1, 1), dots, 2);

  AIM::NearfieldInteraction nf(dots, nullptr, g, 4, 1, 1, grid);

  for(const auto &pair : nf.build_pair_list()) {
    std::cout << pair.first << " " << pair.second << std::endl;
  }


}

BOOST_AUTO_TEST_SUITE_END()  // NEARFIELD
