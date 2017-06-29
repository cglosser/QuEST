#ifndef BLOCH_RHS_H
#define BLOCH_RHS_H

#include <Eigen/Dense>
#include <algorithm>
#include <vector>

#include "../../interactions/interaction.h"
#include "rhs.h"
namespace Integrator {
  class BlochRHS;
}

class Integrator::BlochRHS : public Integrator::RHS<Eigen::Vector2cd> {
 public:
  BlochRHS(const double,
           const std::shared_ptr<History<Eigen::Vector2cd>> &,
           std::vector<std::shared_ptr<Interaction>>,
           std::vector<BlochFunctionType>);
  void evaluate(const int) const override;

 private:
  int num_solutions;
  std::vector<std::shared_ptr<Interaction>> interactions;
  std::vector<BlochFunctionType> rhs_functions;
};

#endif
