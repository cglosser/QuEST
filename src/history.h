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
