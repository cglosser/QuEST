#ifndef HISTORY_H
#define HISTORY_H

#include <boost/multi_array.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>

namespace Integrator {
  template <class soltype>
  class History;

  template <class soltype>
  using soltype_array = boost::multi_array<soltype, 3>;
}

template <class soltype>
class Integrator::History {
 public:
  History(const int, const int, const int);
  soltype_array<soltype> history;

 private:
};

template <class soltype>
Integrator::History<soltype>::History(const int num_particles, const int window,
                                      const int num_timesteps)
    : history(
          boost::extents[num_particles]
                        [typename soltype_array<soltype>::extent_range(-window, num_timesteps)]
                        [2])

{
}


//void History::write_history(
    //const std::shared_ptr<const History::HistoryArray> &history,
    //const std::string &filename, const int n [> = 0 <])
//{
  //std::ofstream output(filename);
  //output << std::setprecision(12) << std::scientific;

  //const int max_t =
      //(n != 0) ? n : history->shape()[1] + history->index_bases()[1];

  //for(int t = 0; t < max_t; ++t) {
    //for(int sol_idx = 0; sol_idx < static_cast<int>(history->shape()[0]);
        //++sol_idx) {
      //output << (*history)[sol_idx][t][0].transpose() << " ";
    //}
    //output << std::endl;
  //}
//}

#endif
