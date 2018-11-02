#ifndef LOGGING_H
#define LOGGING_H

#include <chrono>
#include <fstream>
#include <iomanip>

enum class log_level_t {
  LOG_NOTHING,
  LOG_CRITICAL,
  LOG_ERROR,
  LOG_WARNING,
  LOG_INFO,
  LOG_DEBUG
};

class DurationLogger {
 private:
  std::ofstream log_file;
  std::chrono::high_resolution_clock::time_point t1, t2;

 public:
  DurationLogger(const std::string &fname = "timing.dat")
      : log_file(fname), t1{std::chrono::high_resolution_clock::now()}
  {
  }

  void log_event(const std::string &event)
  {
    t2 = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double, std::micro> duration =
        std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

    log_file << std::scientific << event << " " << duration.count() << "us"
             << std::endl;

    t1 = std::chrono::high_resolution_clock::now();
  }
};

#endif
