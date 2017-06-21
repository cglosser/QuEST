#ifndef HISTORY_H
#define HISTORY_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <memory>

namespace History {
  typedef Eigen::Vector2cd soltype;
  typedef boost::multi_array<soltype, 3> HistoryArray;

  std::shared_ptr<HistoryArray> make_shared_history(const int, const int,
                                                    const int);
  bool isfinite(const soltype &);
  void write_history(const std::shared_ptr<const HistoryArray> &,
                     const std::string &, const int = 0);
}

#endif
#include "history.h"

std::shared_ptr<History::HistoryArray> History::make_shared_history(
    const int num_particles, const int window, const int num_timesteps)
{
  const History::HistoryArray::extent_range time_range(-window, num_timesteps);
  return std::make_shared<History::HistoryArray>(
      boost::extents[num_particles][time_range][2]);
}

bool History::isfinite(const soltype &sol)
{
  return std::isfinite(sol[0].real()) && std::isfinite(sol[0].imag()) &&
         std::isfinite(sol[1].real()) && std::isfinite(sol[1].imag());
}

void History::write_history(
    const std::shared_ptr<const History::HistoryArray> &history,
    const std::string &filename, const int n /* = 0 */)
{
  std::ofstream output(filename);
  output << std::setprecision(12) << std::scientific;

  const int max_t =
      (n != 0) ? n : history->shape()[1] + history->index_bases()[1];

  for(int t = 0; t < max_t; ++t) {
    for(int sol_idx = 0; sol_idx < static_cast<int>(history->shape()[0]);
        ++sol_idx) {
      output << (*history)[sol_idx][t][0].transpose() << " ";
    }
    output << std::endl;
  }
}
