#ifndef LLG_RHS_H
#define LLG_RHS_H

#include "rhs.h"
#include <Eigen/Dense>

namespace Integrator {
  class LLG_RHS;
}

typedef Eigen::Vector3d vec3d;

class Integrator::LLG_RHS : public RHS<vec3d> {
 public:
  LLG_RHS(const double dt,
          const std::shared_ptr<Integrator::History<vec3d>> &history)
      : RHS(dt, history){};
  void evaluate(const int) const;

 private:
};

void Integrator::LLG_RHS::evaluate(const int n) const
{  // constants subject to change
  const double alpha = 1;
  const double gamma0 = 1;
  const double time = n * dt;
  const double gamma = gamma0 / (1 + std::pow(alpha,2)); //might change
//  vec3d hfield = hist_interaction + pulse_interaction; <-- seudo-code
  vec3d hfield(0,0,0);

  for(int i = 0; static_cast<int>(history->array.shape()[0]); ++i) {
    vec3d mxh = history->array[i][n][0].cross(hfield);
    history->array[i][n][1] = -gamma * mxh - gamma * alpha /
	                      history->array[i][n][0].norm() *
                              history->array[i][n][0].cross(mxh);
  };
};
	   
#endif
