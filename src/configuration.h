#ifndef CONFIG_H
#define CONFIG_H

#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>

struct Configuration{
  size_t num_particles;
  double c0, hbar;
  double simulation_time;
};

// This Exception really just allows parse_configs to bail out early if it
// finds a --help or --version flag.
struct CommandLineException : public std::exception {
  const char *what() const throw() {return "Have a nice day! :)";}
};

boost::program_options::variables_map parse_configs(int argc, char *argv[]);
void populate_universe(boost::program_options::variables_map const &);

extern Configuration config;

#endif
