#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

#include "interactions/AIM/aim_interaction.h"
#include "interactions/direct_interaction.h"
#include "math_utils.h"

using dbl = std::numeric_limits<double>;

class Gaussian {
 public:
  Gaussian(double mu, double sigma) : mu(mu), sigma(sigma){};
  double operator()(double t) const { return Math::gaussian((t - mu) / sigma); }
  double deriv(double t) const
  {
    const double arg = (t - mu) / sigma;
    return -arg * Math::gaussian(arg) / sigma;
  }

 private:
  double mu, sigma;
};

double analytic_solution(double c,
                         double t,
                         const Eigen::Vector3d &dr,
                         const Gaussian &f)
{
  return -dr(2) *
         (c * f(t - dr.norm() / c) + dr.norm() * f.deriv(t - dr.norm() / c)) /
         (c * std::pow(dr.norm(), 3));
}

int main()
{
  // Set physical parameters
  const int interpolation_order = 4, expansion_order = 1;

  const int num_steps = 1024;
  const double c = 1, dt = 1, total_time = dt * num_steps;

  const double wavelength = 10 * c * dt;

  double h;
  std::cin >> h;

  const Eigen::Array3d spacing = Eigen::Array3d(1, 1, 1) * c * dt * h;

  auto dots = std::make_shared<DotVector>();

  const int n_pts = 32;
  for(int i = 0; i < n_pts; ++i) {
    dots->push_back(
        QuantumDot(Eigen::Vector3d::Random() * c * dt / 10, {1, 0, 0}));

    dots->push_back(QuantumDot(Eigen::Vector3d::Random() * c * dt / 10 +
                                   Eigen::Vector3d(0, 0, 2 * wavelength),
                               {1, 0, 0}));
  }

  auto grid = AIM::Grid(spacing, expansion_order, *dots);

  // Construct Gaussian history
  const int num_dots = dots->size();
  auto history = std::make_shared<Integrator::History<Eigen::Vector2cd>>(
      num_dots, 10, num_steps);

  const Gaussian source(total_time / 2.0, total_time / 12.0);
  for(int i = 0; i < n_pts; ++i) {
    for(int t = -10; t < num_steps; ++t) {
      history->array_[i][t][0](RHO_01) = source(t * dt);
    }
  }

  // Construct interactions
  Propagation::Laplace<cmplx> prop;
  int transit = grid.max_transit_steps(c, dt) + interpolation_order;

  AIM::Interaction aim(dots, history, prop, spacing, interpolation_order,
                       expansion_order, 1, c, dt, AIM::Normalization::Laplace(),
                       Projector::Derivative<cmplx>(transit));
  DirectInteraction direct(dots, history, prop, interpolation_order, c, dt);

  std::ofstream fd("output.dat");
  fd.precision(17);

  std::cout << "Shape: " << grid.shape().transpose() << std::endl;

  for(int t = 0; t < num_steps; ++t) {
    if(t % 100 == 0) std::cout << t << std::endl;

    Eigen::RowVectorXd analytic = Eigen::VectorXd::Zero(n_pts);

    for(int src = 0; src < n_pts; ++src) {
      for(int obs = 0; obs < n_pts; ++obs) {
        analytic(obs) += analytic_solution(
            c, t * dt,
            dots->at(n_pts + obs).position() - dots->at(src).position(),
            source);
      }
    }

    auto aux = aim.evaluate(t);

    Eigen::Map<Eigen::RowVectorXcd> b(aux.data() + n_pts, n_pts);
    // fd << (analytic - b).norm() / analytic.norm() << std::endl;

    fd << b.real() << " " << analytic << std::endl;
  }

  return 0;
}
