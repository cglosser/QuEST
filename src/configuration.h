#ifndef CONFIG_H
#define CONFIG_H

#include <boost/program_options.hpp>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

struct Configuration {
  int num_particles;
  int num_timesteps;
  double c0, hbar, mu0, omega;
  double dt, total_time;
  int interpolation_order;
  std::string qd_path, pulse_path;
};

// This Exception really just allows parse_configs to return control
// to main() in an atypical manner (--help or --version flag, config
// file not found, etc.)
struct CommandLineException : public std::exception {
  const char *what() const throw() { return "Have a nice day! :)"; }
};

boost::program_options::variables_map parse_configs(int argc, char *argv[]);

extern Configuration config;

#endif
