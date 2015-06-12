#include <iostream>
#include <fstream>
#include <iterator>
using namespace std;

#include <boost/program_options.hpp>
namespace po = boost::program_options;

// A helper function to simplify the main part.
template<class T>
ostream& operator<<(ostream& os, const vector<T>& v) {
  copy(v.begin(), v.end(), ostream_iterator<T>(os, " "));
  return os;
}

int main(int argc, char *argv[]) {
  try {
    double c0, hbar;
    string config_file;

    // Declare a group of options that will be
    // allowed only on command line
    po::options_description generic("Command line options");
    generic.add_options()
      ("help", "print this help message")
      ("version,v", "print version string")
      ("config,c", po::value<string>(&config_file)->default_value("input.cfg"),
        "path to configuration file")
      ;

    // Declare a group of options that will be
    // allowed both on command line and in
    // config file
    po::options_description config("Physical constants");
    config.add_options()
      ("constants.c0", po::value<double>(&c0),
        "speed of light in vacuum")
      ("constants.hbar", po::value<double>(&hbar),
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
      return 0;
    }

    if (vm.count("version")) {
      cout << "ElectroMagnetics Research Gadget, version 0" << endl;
      return 0;
    }

    ifstream ifs(config_file.c_str());
    if (!ifs) {
      cerr << "ERROR: " << config_file << " not found" << endl;
      return 0;
    } else {
      po::store(po::parse_config_file(ifs, config_file_options), vm);
      po::notify(vm);
    }

    cout << "speed of light: " << c0 << endl;
    cout << "          hbar: " << hbar << endl;
  } catch(exception &e) {
    cerr << e.what() << "\n";
    return 1;
  }
  return 0;
}
