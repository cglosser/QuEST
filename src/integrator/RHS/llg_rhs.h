#ifndef LLG_RHS_H
#define LLG_RHS_H

#include <Eigen/Dense>

#include "../../interactions/interaction.h"
#include "rhs.h"

namespace Integrator {
  class LLG_RHS;
}

class Integrator::LLG_RHS : public RHS<Eigen::Vector3d> {
 public:
  LLG_RHS(const double,
          const std::shared_ptr<Integrator::History<soltype>> &,
          std::vector<std::shared_ptr<Interaction>>,
          rhs_func_vector);
  void evaluate(const int) const override;

 private:
  int num_solutions;
  std::vector<std::shared_ptr<Interaction>> interactions;
  rhs_func_vector rhs_functions;
};

#endif
