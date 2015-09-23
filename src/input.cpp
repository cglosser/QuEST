#include "input.h"

using namespace std;
namespace po = boost::program_options;

po::variables_map parse_configs(int argc, char *argv[]) {
  string config_path, domain_keyword;

  po::options_description cmd_line_desc("Command line options");
  cmd_line_desc.add_options()
    ("help", "print this help message")
    ("version,v", "print version string")
    ("config,c", po::value<string>(&config_path)->default_value("input.cfg"), "path to configuration file")
  ;

  po::options_description system_parameters("System parameters");
  system_parameters.add_options()
    ("parameters.dimensions",      po::value<size_t>(&Universe.dimensions)->default_value(3), "number of spatial dimensions")
    ("parameters.num_particles",   po::value<size_t>(&Universe.num_particles)->required(), "number of particles in the system")
    ("parameters.simulation_time", po::value<double>(&Universe.simulation_time)->required(), "total (time-domain) simulation duration")
  ;

  po::options_description constants("Physical constants");
  constants.add_options()
    ("constants.c0",   po::value<double>(&Universe.c0)->default_value(1.0), "speed of light in vacuum")
    ("constants.hbar", po::value<double>(&Universe.hbar)->default_value(1.0), "reduced Planck constant")
  ;

  po::options_description cmdline_options;
  cmdline_options.add(cmd_line_desc);

  po::options_description config_file_options;
  config_file_options.add(system_parameters);
  config_file_options.add(constants);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    po::options_description visible("EMRG options");
    visible.add(cmd_line_desc).add(system_parameters).add(constants);
    cout << visible << "\n";
    throw CommandLineException();
  }

  if (vm.count("version")) {
    cout << "ElectroMagnetics Research Gadget, version 0" << endl;
    throw CommandLineException();
  }

  ifstream ifs(config_path.c_str());
  if (!ifs) {
    cerr << "ERROR: " << config_path << " not found" << endl;
  } else {
    po::store(po::parse_config_file(ifs, config_file_options), vm);
    po::notify(vm);
  }

  return vm;
}
