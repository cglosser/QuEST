#ifndef NORMALIZATION_H
#define NORMALIZATION_H

namespace AIM {
  namespace Normalization {
    using SpatialNorm = std::function<cmplx(const Eigen::Vector3d &)>;
    const SpatialNorm unit = [](__attribute__((unused))
                                const Eigen::Vector3d &v) { return 1; };

    class Laplace {
     public:
      Laplace(const double alpha = 1) : alpha(alpha){};
      double operator()(const Eigen::Vector3d &dr) const
      {
        return alpha / dr.norm();
      }

     private:
      double alpha;
    };

    class Helmholtz {
     public:
      Helmholtz(const double k = 0, const double alpha = 1)
          : k(k), alpha(alpha){};
      cmplx operator()(const Eigen::Vector3d &dr) const
      {
        const double R = dr.norm();
        return std::exp(-iu * k * R) * alpha / R;
      }

     private:
      double k, alpha;
    };
  }
}

#endif
