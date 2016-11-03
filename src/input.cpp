#include "input.h"

using namespace std;
namespace po = boost::program_options;

po::variables_map parse_configs(int argc, char *argv[]) {
  string config_path, domain_keyword;

  po::options_description cmd_line_description("Command line options");
  cmd_line_description.add_options()
    ("help", "print this help message")
    ("version,v", "print version string")
    ("config,c", po::value<string>(&config_path)->default_value("input.cfg"), "path to configuration file")
  ;

  po::options_description file_description("System parameters");
  file_description.add_options()
    ("parameters.num_particles",   po::value<size_t>(&Universe.num_particles)->required(), "number of particles in the system")
    ("parameters.simulation_time", po::value<double>(&Universe.simulation_time)->required(), "total (time-domain) simulation duration")

    ("constants.c0",   po::value<double>(&Universe.c0)->default_value(1.0), "speed of light in vacuum")
    ("constants.hbar", po::value<double>(&Universe.hbar)->default_value(1.0), "reduced Planck constant")
  ;

  po::options_description cmdline_options, file_options;
  cmdline_options.add(cmd_line_description);
  file_options.add(file_description);

  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(cmdline_options).run(), vm);
  po::notify(vm);

  if (vm.count("help")) {
    po::options_description visible("EMRG options");
    visible.add(cmd_line_description).add(file_description);
    cout << visible << "\n";
    throw CommandLineException();
  }

  if (vm.count("version")) {
    cout << "ElectroMagnetics Research Gadget, version 0" << endl;
    cout << "Compiled with " << __VERSION__ << " on " << __TIMESTAMP__ << endl;
    throw CommandLineException();
  }

  ifstream ifs(config_path.c_str());
  if (!ifs) {
    cerr << "ERROR: " << config_path << " not found" << endl;
  } else {
    po::store(po::parse_config_file(ifs, file_options), vm);
    po::notify(vm);
  }

  return vm;
}
