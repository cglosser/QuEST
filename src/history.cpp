#include "history.h"

std::shared_ptr<History::HistoryArray> History::make_shared_history(
    const int num_particles, const int window, const int num_timesteps)
{
  const History::HistoryArray::extent_range time_range(-window, num_timesteps);
  return std::make_shared<History::HistoryArray>(
      boost::extents[num_particles][time_range][2]);
}
