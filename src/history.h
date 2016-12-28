#ifndef HISTORY_H
#define HISTORY_H

#include <Eigen/Dense>
#include <boost/multi_array.hpp>

namespace History {
  typedef Eigen::Vector2cd soltype;
  typedef boost::multi_array<soltype, 3> HistoryArray;
}

#endif
