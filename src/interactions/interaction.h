#ifndef INTERACTION_H
#define INTERACTION_H

#include <Eigen/Dense>
#include <memory>

#include "../quantum_dot.h"

class Interaction {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  Interaction(const std::shared_ptr<const DotVector> dots, const double dt)
      : dots(std::move(dots)), results(dots ? dots->size() : 0), dt(dt){};
  const cmplx &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;
  virtual ~Interaction() {};

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
  double dt;
};

#endif
