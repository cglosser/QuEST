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
