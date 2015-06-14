#include "input.h"

using namespace std;
namespace po = boost::program_options;

po::variables_map parse_configs(int argc, char *argv[]) {
  string config_file;

  //Command-line only options; really just for printing program info
  po::options_description generic("Command line options");
  generic.add_options()
    ("help", "print this help message")
    ("version,v", "print version string")
    ("config,c", po::value<string>(&config_file)->default_value("input.cfg"),
      "path to configuration file")
    ;

  //Physical constants to read from a config file
  po::options_description config("Physical constants");
  config.add_options()
    ("constants.c0", po::value<double>()->default_value(1.0),
      "speed of light in vacuum")
    ("constants.hbar", po::value<double>()->default_value(1.0),
      "reduced Planck constant")
    ;

  po::options_description cmdline_options;
  cmdline_options.add(generic);

  po::options_description config_file_options;
  config_file_options.add(config);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    po::options_description visible("EMRG options");
    visible.add(generic).add(config);
    cout << visible << "\n";
    throw SilentException();
  }

  if (vm.count("version")) {
    cout << "ElectroMagnetics Research Gadget, version 0" << endl;
    throw SilentException();
  }

  ifstream ifs(config_file.c_str());
  if (!ifs) {
    cerr << "ERROR: " << config_file << " not found" << endl;
  } else {
    po::store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);
  }

  return vm;
}

void populate_universe(po::variables_map const &vm) {
  Universe.c0   = vm["constants.c0"].as<double>();
  Universe.hbar = vm["constants.hbar"].as<double>();
}
