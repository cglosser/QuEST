#include "../src/interactions/history_interaction.h"
#include <boost/test/unit_test.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(history_interaction)

BOOST_AUTO_TEST_CASE(greens_function_coefficients)
{
  Eigen::Vector3d pos1(0, 0, 0);
  Eigen::Vector3d pos2(1, 1, 1);
  const std::pair<double, double> damping(1, 1);
  const double freq = 2278.9013;
  const Eigen::Vector3d dip(1, 2, 3);

  QuantumDot qd1(pos1, freq, damping, dip);
  QuantumDot qd2(pos2, freq, damping, dip);

  DotVector dots_vec = {qd1,qd2};
  auto dots = std::make_shared<DotVector>(dots_vec);

  auto history(History::make_shared_history(2, 22, 1));
  for(int dot_idx=0; dot_idx<2; ++dot_idx) {
    for(int time_idx=-22; time_idx<=0; ++time_idx) {
      (*history)[dot_idx][time_idx][0] = History::soltype(1, 0);
    }
  }

  auto dyadic(std::make_shared<GreenFunction::Dyadic>(
	      GreenFunction::Dyadic(1, 2, 1)));
  auto hist_inter = HistoryInteraction(dots, history, dyadic, 1, 1, 1);

  auto result = hist_inter.evaluate(1);
  std::cout << result << std::endl;
}

BOOST_AUTO_TEST_SUITE_END()
