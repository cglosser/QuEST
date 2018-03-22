#pragma once

#include <Eigen/Dense>
#include <boost/program_options.hpp>
#include <cmath>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

struct Configuration {
  // Parameters
  int num_particles, num_timesteps, interpolation_order;
  double dt, total_time;
  enum class SIMULATION_TYPE { DIRECT, FAST } sim_type;

  // File paths
  std::string qd_path, pulse_path;

  // Constants
  double c0, mu0, hbar, laser_frequency;

  // AIM
  int expansion_order, border;
  Eigen::Vector3d grid_spacing;
};

// This Exception really just allows parse_configs to return control
// to main() in an atypical manner (--help or --version flag, config
// file not found, etc.)
struct CommandLineException : public std::exception {
  const char *what() const throw() { return "Have a nice day! :)"; }
};

boost::program_options::variables_map parse_configs(int argc, char *argv[]);

extern Configuration config;
