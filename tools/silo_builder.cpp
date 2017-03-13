#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <silo.h>
#include <sstream>
#include <vector>

using std::cout;
using std::endl;

constexpr int nmax = 1024;

typedef Eigen::Matrix<double, nmax, 3> PosArray; // Column major!

class SiloFile {
public:
  SiloFile(const std::string &fname) { db_file_ptr = DBCreate(fname.data(),
      DB_CLOBBER, DB_LOCAL, "bloch", DB_HDF5); }

  ~SiloFile() { DBClose(db_file_ptr); }

  void write_point_mesh(const PosArray &coords, const std::string &mesh_id) {
    double data_ptr[coords.size()]; Eigen::Map<PosArray>(data_ptr,
        coords.rows(), coords.cols()) = coords;

    double *coord_ptrs[] = {&data_ptr[0], &data_ptr[coords.rows()],
      &data_ptr[2*coords.rows()]};

    DBPutPointmesh(db_file_ptr, mesh_id.data(), 3, coord_ptrs, coords.rows(),
                   DB_DOUBLE, NULL);
  }

  void write_point_data(std::vector<double> &data, const std::string
      &name, const std::string &mesh_name) {

    DBPutPointvar1(db_file_ptr, name.data(), mesh_name.data(), data.data(),
        data.size(), DB_DOUBLE, NULL);
  }

private:
  DBfile *db_file_ptr;
};


PosArray read_coords(const std::string &);

int main() {

  auto xyz(read_coords("/home/connor/Scratch/neighborhood_00/dots_neighborhood00.cfg"));

  SiloFile sf("potato.silo");

  sf.write_point_mesh(xyz, "coordinates");

  SiloFile sf00("potato_00.silo");
  std::vector<double> d(1024, 1);
  sf00.write_point_mesh(xyz, "coordinates");
  sf00.write_point_data(d, "pop", "coordinates");

  SiloFile sf01("potato_01.silo");
  for(auto &i : d) i = 2;
  sf01.write_point_mesh(xyz, "coordinates");
  sf01.write_point_data(d, "pop", "coordinates");

  SiloFile sf02("potato_02.silo");
  for(auto &i : d) i = 3;
  sf02.write_point_mesh(xyz, "coordinates");
  sf02.write_point_data(d, "pop", "coordinates");

  SiloFile sf03("potato_03.silo");
  for(auto &i : d) i = 4;
  sf03.write_point_mesh(xyz, "coordinates");
  sf03.write_point_data(d, "pop", "coordinates");

  SiloFile sf04("potato_04.silo");
  for(auto &i : d) i = 5;
  sf04.write_point_mesh(xyz, "coordinates");
  sf04.write_point_data(d, "pop", "coordinates");

  SiloFile sf05("potato_05.silo");
  for(auto &i : d) i = 6;
  sf05.write_point_mesh(xyz, "coordinates");
  sf05.write_point_data(d, "pop", "coordinates");

  SiloFile sf06("potato_06.silo");
  for(auto &i : d) i = 7;
  sf06.write_point_mesh(xyz, "coordinates");
  sf06.write_point_data(d, "pop", "coordinates");

  SiloFile sf07("potato_07.silo");
  for(auto &i : d) i = 8;
  sf07.write_point_mesh(xyz, "coordinates");
  sf07.write_point_data(d, "pop", "coordinates");

  SiloFile sf08("potato_08.silo");
  for(auto &i : d) i = 9;
  sf08.write_point_mesh(xyz, "coordinates");
  sf08.write_point_data(d, "pop", "coordinates");

  SiloFile sf09("potato_09.silo");
  for(auto &i : d) i = 10;
  sf09.write_point_mesh(xyz, "coordinates");
  sf09.write_point_data(d, "pop", "coordinates");




  return 0;
}

PosArray read_coords(const std::string &fname) {
  std::ifstream ifs(fname);
  PosArray result;
  std::string line;

  for(int i = 0; i < nmax; ++i) {
    std::getline(ifs, line);
    std::istringstream ss(line);
    ss >> result(i,0) >> result(i,1) >> result(i,2);
  }

  return result;
}
