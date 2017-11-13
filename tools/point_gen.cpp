#include <Eigen/Dense>
#include <algorithm>
#include <array>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <vector>

using std::cout;
using std::endl;

typedef std::array<double, 9> Qdot;

std::ostream &operator<<(std::ostream &out, const Qdot &qd)
{
  std::copy(qd.begin(), qd.end(), std::ostream_iterator<double>(out, " "));
  return out;
}

double dist(const Qdot &qd1, const Qdot &qd2)
{
  Eigen::Vector3d dr(qd2[0] - qd1[0], qd2[1] - qd1[1], qd2[2] - qd1[2]);
  return dr.norm();
}

double min_dist(const std::vector<Qdot> &dots)
{
  int idx1 = -1, idx2 = -1;
  double min = 1e100;
  for(size_t i = 0; i < dots.size() - 1; ++i) {
    for(size_t j = i + 1; j < dots.size(); ++j) {
      double dij = dist(dots[i], dots[j]);
      if(dij < min) {
        min = dij;
        idx1 = i;
        idx2 = j;
      }
    }
  }

  cout << "Min dist: " << min << " (" << idx1 << "," << idx2 << ")" << endl;

  return min;
}

int main()
{
  const int num_dots = 10000;
  const double xlen = 2.4, ylen = 2.4, zlen = 0.1;

  const unsigned seed1 =
      std::chrono::system_clock::now().time_since_epoch().count();

  std::default_random_engine generator(seed1);
  std::uniform_real_distribution<double> xdist(-xlen / 2, xlen / 2),
      ydist(-ylen / 2, ylen / 2), zdist(-zlen / 2, zlen / 2);

  std::vector<Qdot> dots;
  dots.reserve(num_dots);

  for(int i = 0; i < num_dots; ++i) {
    Qdot qd = {{xdist(generator), ydist(generator), zdist(generator), 2278.9013,
                10.0, 20.0, 5.2917721e-4, 0.0, 0.0}};
    dots.push_back(qd);
  }

  std::sort(dots.begin(), dots.end(), [](const Qdot &a, const Qdot &b) {
    return a[0] * a[0] + a[1] * a[1] + a[2] * a[2] <
           b[0] * b[0] + b[1] * b[1] + b[2] * b[2];
  });

  min_dist(dots);

  std::ofstream fd("dots.dat");
  fd << std::setprecision(12);

  for(const auto &d : dots) {
    fd << d << endl;
  }

  fd.close();

  return 0;
}
