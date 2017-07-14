#include <boost/test/unit_test.hpp>
#include <cmath>

#include "../src/interactions/aim_interaction.h"

BOOST_AUTO_TEST_SUITE(AIM)

BOOST_AUTO_TEST_SUITE(GridSort)

BOOST_AUTO_TEST_CASE(PointMapping)
{
  const Eigen::Vector3d dip(0, 0, 1);
  const auto damping = std::make_pair(10.0, 20.0);
  Eigen::Vector3d grid_spacing(1, 1, 1);
  auto dots = std::make_shared<DotVector>(DotVector(
      {QuantumDot(Eigen::Vector3d(-0.5, -0.5, 0), 1.0, damping, dip),
       QuantumDot(Eigen::Vector3d(-0.2, 1.3, 0), 1.0, damping, dip),
       QuantumDot(Eigen::Vector3d(0.2, 1.3, 0), 1.0, damping, dip),
       QuantumDot(Eigen::Vector3d(0.4, -0.7, 0), 1.0, damping, dip)}));

  AIM::Grid grid(grid_spacing, dots);

  for(size_t box_idx = 0; box_idx < grid.point_mapping.size(); ++box_idx) {
    std::cout << box_idx << ": ";
    for(size_t particle_idx = 0; particle_idx < grid.point_mapping.at(box_idx).size(); ++particle_idx) {
      std::cout << grid.point_mapping.at(box_idx).at(particle_idx) << " ";
    }
    std::cout << std::endl;
  }

}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE_END()
