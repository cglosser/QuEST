#ifndef INPUT_H
#define INPUT_H

#include <boost/program_options.hpp>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include "universe.h"

boost::program_options::variables_map parse_configs(int argc, char *argv[]);
void populate_universe(boost::program_options::variables_map const &);

#endif
