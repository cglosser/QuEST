#ifndef INPUT_H
#define INPUT_H

#include <boost/program_options.hpp>
#include <exception>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include "universe.h"

// This Exception really just allows parse_configs to bail out early if it
// finds a --help or --version flag.
struct CommandLineException : public std::exception {
  const char *what() const throw() {return "Have a nice day! :)";}
};

boost::program_options::variables_map parse_configs(int argc, char *argv[]);
void populate_universe(boost::program_options::variables_map const &);

#endif
