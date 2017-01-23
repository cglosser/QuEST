#ifndef INTERACTION_H
#define INTERACTION_H

#include <memory>
#include <Eigen/Dense>

#include "../quantum_dot.h"

class Interaction {
 public:
  typedef Eigen::Array<cmplx, Eigen::Dynamic, 1> ResultArray;

  Interaction(const std::shared_ptr<const DotVector> &dots)
      : dots(dots), results(dots->size()){};
  const cmplx &operator[](const int i) const { return results[i]; }
  virtual const ResultArray &evaluate(const int) = 0;

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultArray results;
};

#endif
