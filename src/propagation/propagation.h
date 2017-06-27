#ifndef PROPAGATION_H
#define PROPAGATION_H

#include <Eigen/Dense>
#include <Eigen/StdVector>
#include <vector>

#include "../quantum_dot.h"

namespace Propagation {
  template <class fieldtype>
  class FieldEvaluator;
}

template <class scalarFieldType>
class Propagation::FieldEvaluator {
 public:
  typedef Eigen::Matrix<scalarFieldType, 1, 3> fieldtype;
  typedef std::vector<fieldtype, Eigen::aligned_allocator<fieldtype>>
      ResultTable;

  FieldEvaluator(const std::shared_ptr<const DotVector> &dots)
      : dots(dots), results(dots->size()){};

  const fieldtype &operator[](const int i) const { return results[i]; }
  virtual const ResultTable &evaluate(const int) = 0;

 protected:
  std::shared_ptr<const DotVector> dots;
  ResultTable results;
};

#endif
