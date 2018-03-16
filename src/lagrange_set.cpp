#include "lagrange_set.h"

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const int order)
    : evaluations(boost::extents[Interpolation::NUM_DERIVATIVES][order + 1]),
      order_(order)
{
}

Interpolation::UniformLagrangeSet::UniformLagrangeSet(const double x,
                                                      const int order,
                                                      const double dt /* = 1 */)
    : UniformLagrangeSet(order)
{
  assert(x > 0);  // Don't extrapolate!
  evaluate_derivative_table_at_x(x, dt);
}

void Interpolation::UniformLagrangeSet::evaluate_derivative_table_at_x(
    const double x, const double dt /* = 1 */)
{
  for(int basis_id = 0; basis_id <= order_; ++basis_id) {
    double d0_product = 1, d1_sum = 0, d2_sum = 0;
    for(int m = 0; m <= order_; ++m) {
      if(m == basis_id) continue;

      double numer = (x - m);

      d0_product *= numer / (basis_id - m);

      if(numer != 0) {
        d1_sum -= 1 / numer;  // Note the minus sign!
        d2_sum -= std::pow(numer, -2);
      }
    }

    evaluations[0][basis_id] = d0_product;
    evaluations[1][basis_id] = (evaluations[0][basis_id] * d1_sum);
    evaluations[2][basis_id] = (evaluations[0][basis_id] * d2_sum +
         std::pow(evaluations[1][basis_id], 2) / evaluations[0][basis_id]);
  }

  for(int i = 0; i <= order_; ++i) {
    evaluations[1][i] *= std::pow(dt, -1);
    evaluations[2][i] *= std::pow(dt, -2);
  }
}

//clang-format off
double PolyTree(const int order, const int derivative, const double x)
{
  const int id = static_cast<int>(std::ceil(x));
  const double x_sq = x * x;

  switch(order) {
    case(0):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1;
            default: return 0;
          }

        default: return 0;
      }

    case(1):
      switch(derivative) {
        case(0):
          switch(id) {
            case 0: return 1 + x;
            case 1: return 1 - x;
            default: return 0;
          }

        case(1):
          switch(id) {
            case 0: return 1;
            case 1: return -1;
            default: return 0;
          }

        default: return 0;
      }

    case(2):
      switch(derivative) {
        case(0):
          switch(id) {
            case 0: return 1 + (1.5 + x / 2.) * x;
            case 1: return 1 - x_sq;
            case 2: return 1 + (-1.5 + x / 2.) * x;
            default: return 0;
          }

        case(1):
          switch(id) {
            case 0: return 1.5 + x;
            case 1: return -2 * x;
            case 2: return -1.5 + x;
            default: return 0;
          }

        case(2):
          switch(id) {
            case 0: return 1;
            case 1: return -2;
            case 2: return 1;
            default: return 0;
          }

        default: return 0;
      }

    case(3):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (1.8333333333333333 + (1 + x / 6.) * x);
            case 1: return 1 + x * (0.5 + (-1 - x / 2.) * x);
            case 2: return 1 + x * (-0.5 + (-1 + x / 2.) * x);
            case 3: return 1 + x * (-1.8333333333333333 + (1 - x / 6.) * x);
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 1.8333333333333333 + (2 + x / 2.) * x;
            case 1: return 0.5 + (-2 - (3 * x) / 2.) * x;
            case 2: return -0.5 + x * (-2 + (3 * x) / 2.);
            case 3: return -1.8333333333333333 + (2 - x / 2.) * x;
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 2 + x;
            case 1: return -2 - 3 * x;
            case 2: return -2 + 3 * x;
            case 3: return 2 - x;
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 1;
            case 1: return -3;
            case 2: return 3;
            case 3: return -1;
            default: return 0;
          }

        default: return 0;
      }

    case(4):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (2.0833333333333335 + x * (1.4583333333333333 + (0.4166666666666667 + x / 24.) * x));
            case 1: return 1 + x * (0.8333333333333334 + x * (-0.8333333333333334 + (-0.8333333333333334 - x / 6.) * x));
            case 2: return 1 + x_sq * (-1.25 + x_sq / 4.);
            case 3: return 1 + x * (-0.8333333333333334 + x * (-0.8333333333333334 + (0.8333333333333334 - x / 6.) * x));
            case 4: return 1 + x * (-2.0833333333333335 + x * (1.4583333333333333 + (-0.4166666666666667 + x / 24.) * x));
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 2.0833333333333335 + x * (2.9166666666666665 + (1.25 + x / 6.) * x);
            case 1: return 0.8333333333333334 + x * (-1.6666666666666667 + (-2.5 - (2 * x) / 3.) * x);
            case 2: return x * (-2.5 + x_sq);
            case 3: return -0.8333333333333334 + x * (-1.6666666666666667 + (2.5 - (2 * x) / 3.) * x);
            case 4: return -2.0833333333333335 + x * (2.9166666666666665 + (-1.25 + x / 6.) * x);
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 2.9166666666666665 + (2.5 + x / 2.) * x;
            case 1: return -1.6666666666666667 + (-5 - 2 * x) * x;
            case 2: return -2.5 + 3 * x_sq;
            case 3: return -1.6666666666666667 + (5 - 2 * x) * x;
            case 4: return 2.9166666666666665 + (-2.5 + x / 2.) * x;
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 2.5 + x;
            case 1: return -5 - 4 * x;
            case 2: return 6 * x;
            case 3: return 5 - 4 * x;
            case 4: return -2.5 + x;
            default: return 0;
          }

        case 4:
          switch(id) {
            case 0: return 1;
            case 1: return -4;
            case 2: return 6;
            case 3: return -4;
            case 4: return 1;
            default: return 0;
          }

        default: return 0;
      }

    case(5):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (2.283333333333333 + x * (1.875 + x * (0.7083333333333334 + (0.125 + x / 120.) * x)));
            case 1: return 1 + x * (1.0833333333333333 + x * (-0.625 + x * (-1.0416666666666667 + (-0.375 - x / 24.) * x)));
            case 2: return 1 + x * (0.3333333333333333 + x * (-1.25 + x * (-0.4166666666666667 + (0.25 + x / 12.) * x)));
            case 3: return 1 + x * (-0.3333333333333333 + x * (-1.25 + x * (0.4166666666666667 + (0.25 - x / 12.) * x)));
            case 4: return 1 + x * (-1.0833333333333333 + x * (-0.625 + x * (1.0416666666666667 + (-0.375 + x / 24.) * x)));
            case 5: return 1 + x * (-2.283333333333333 + x * (1.875 + x * (-0.7083333333333334 + (0.125 - x / 120.) * x)));
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 2.283333333333333 + x * (3.75 + x * (2.125 + (0.5 + x / 24.) * x));
            case 1: return 1.0833333333333333 + x * (-1.25 + x * (-3.125 + (-1.5 - (5 * x) / 24.) * x));
            case 2: return 0.3333333333333333 + x * (-2.5 + x * (-1.25 + (1 + (5 * x) / 12.) * x));
            case 3: return -0.3333333333333333 + x * (-2.5 + x * (1.25 + (1 - (5 * x) / 12.) * x));
            case 4: return -1.0833333333333333 + x * (-1.25 + x * (3.125 + (-1.5 + (5 * x) / 24.) * x));
            case 5: return -2.283333333333333 + x * (3.75 + x * (-2.125 + (0.5 - x / 24.) * x));
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 3.75 + x * (4.25 + (1.5 + x / 6.) * x);
            case 1: return -1.25 + x * (-6.25 + (-4.5 - (5 * x) / 6.) * x);
            case 2: return -2.5 + x * (-2.5 + x * (3 + (5 * x) / 3.));
            case 3: return -2.5 + x * (2.5 + (3 - (5 * x) / 3.) * x);
            case 4: return -1.25 + x * (6.25 + (-4.5 + (5 * x) / 6.) * x);
            case 5: return 3.75 + x * (-4.25 + (1.5 - x / 6.) * x);
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 4.25 + (3 + x / 2.) * x;
            case 1: return -6.25 + (-9 - (5 * x) / 2.) * x;
            case 2: return -2.5 + x * (6 + 5 * x);
            case 3: return 2.5 + (6 - 5 * x) * x;
            case 4: return 6.25 + x * (-9 + (5 * x) / 2.);
            case 5: return -4.25 + (3 - x / 2.) * x;
            default: return 0;
          }

        case 4:
          switch(id) {
            case 0: return 3 + x;
            case 1: return -9 - 5 * x;
            case 2: return 6 + 10 * x;
            case 3: return 6 - 10 * x;
            case 4: return -9 + 5 * x;
            case 5: return 3 - x;
            default: return 0;
          }

        case 5:
          switch(id) {
            case 0: return 1;
            case 1: return -5;
            case 2: return 10;
            case 3: return -10;
            case 4: return 5;
            case 5: return -1;
            default: return 0;
          }

        default: return 0;
      }

    case(6):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (2.45 + x * (2.2555555555555555 + x * (1.0208333333333333 + x * (0.24305555555555555 + (0.029166666666666667 + x / 720.) * x))));
            case 1: return 1 + x * (1.2833333333333334 + x * (-0.4083333333333333 + x * (-1.1666666666666667 + x * (-0.5833333333333334 + (-0.11666666666666667 - x / 120.) * x))));
            case 2: return 1 + x * (0.5833333333333334 + x * (-1.1666666666666667 + x * (-0.7291666666666666 + x * (0.14583333333333334 + (0.14583333333333334 + x / 48.) * x))));
            case 3: return 1 + x_sq * (-1.3611111111111112 + x_sq * (0.3888888888888889 - x_sq / 36.));
            case 4: return 1 + x * (-0.5833333333333334 + x * (-1.1666666666666667 + x * (0.7291666666666666 + x * (0.14583333333333334 + (-0.14583333333333334 + x / 48.) * x))));
            case 5: return 1 + x * (-1.2833333333333334 + x * (-0.4083333333333333 + x * (1.1666666666666667 + x * (-0.5833333333333334 + (0.11666666666666667 - x / 120.) * x))));
            case 6: return 1 + x * (-2.45 + x * (2.2555555555555555 + x * (-1.0208333333333333 + x * (0.24305555555555555 + (-0.029166666666666667 + x / 720.) * x))));
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 2.45 + x * (4.511111111111111 + x * (3.0625 + x * (0.9722222222222222 + (0.14583333333333334 + x / 120.) * x)));
            case 1: return 1.2833333333333334 + x * (-0.8166666666666667 + x * (-3.5 + x * (-2.3333333333333335 + (-0.5833333333333334 - x / 20.) * x)));
            case 2: return 0.5833333333333334 + x * (-2.3333333333333335 + x * (-2.1875 + x * (0.5833333333333334 + (0.7291666666666666 + x / 8.) * x)));
            case 3: return x * (-2.7222222222222223 + x_sq * (1.5555555555555556 - x_sq / 6.));
            case 4: return -0.5833333333333334 + x * (-2.3333333333333335 + x * (2.1875 + x * (0.5833333333333334 + (-0.7291666666666666 + x / 8.) * x)));
            case 5: return -1.2833333333333334 + x * (-0.8166666666666667 + x * (3.5 + x * (-2.3333333333333335 + (0.5833333333333334 - x / 20.) * x)));
            case 6: return -2.45 + x * (4.511111111111111 + x * (-3.0625 + x * (0.9722222222222222 + (-0.14583333333333334 + x / 120.) * x)));
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 4.511111111111111 + x * (6.125 + x * (2.9166666666666665 + (0.5833333333333334 + x / 24.) * x));
            case 1: return -0.8166666666666667 + x * (-7 + x * (-7 + (-2.3333333333333335 - x / 4.) * x));
            case 2: return -2.3333333333333335 + x * (-4.375 + x * (1.75 + (2.9166666666666665 + (5 * x) / 8.) * x));
            case 3: return -2.7222222222222223 + x_sq * (4.666666666666667 - (5 * x_sq) / 6.);
            case 4: return -2.3333333333333335 + x * (4.375 + x * (1.75 + (-2.9166666666666665 + (5 * x) / 8.) * x));
            case 5: return -0.8166666666666667 + x * (7 + x * (-7 + (2.3333333333333335 - x / 4.) * x));
            case 6: return 4.511111111111111 + x * (-6.125 + x * (2.9166666666666665 + (-0.5833333333333334 + x / 24.) * x));
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 6.125 + x * (5.833333333333333 + (1.75 + x / 6.) * x);
            case 1: return -7 + x * (-14 + (-7 - x) * x);
            case 2: return -4.375 + x * (3.5 + x * (8.75 + (5 * x) / 2.));
            case 3: return x * (9.333333333333334 - (10 * x_sq) / 3.);
            case 4: return 4.375 + x * (3.5 + x * (-8.75 + (5 * x) / 2.));
            case 5: return 7 + x * (-14 + (7 - x) * x);
            case 6: return -6.125 + x * (5.833333333333333 + (-1.75 + x / 6.) * x);
            default: return 0;
          }

        case 4:
          switch(id) {
            case 0: return 5.833333333333333 + (3.5 + x / 2.) * x;
            case 1: return -14 + (-14 - 3 * x) * x;
            case 2: return 3.5 + x * (17.5 + (15 * x) / 2.);
            case 3: return 9.333333333333334 - 10 * x_sq;
            case 4: return 3.5 + x * (-17.5 + (15 * x) / 2.);
            case 5: return -14 + (14 - 3 * x) * x;
            case 6: return 5.833333333333333 + (-3.5 + x / 2.) * x;
            default: return 0;
          }

        case 5:
          switch(id) {
            case 0: return 3.5 + x;
            case 1: return -14 - 6 * x;
            case 2: return 17.5 + 15 * x;
            case 3: return -20 * x;
            case 4: return -17.5 + 15 * x;
            case 5: return 14 - 6 * x;
            case 6: return -3.5 + x;
            default: return 0;
          }

        case 6:
          switch(id) {
            case 0: return 1;
            case 1: return -6;
            case 2: return 15;
            case 3: return -20;
            case 4: return 15;
            case 5: return -6;
            case 6: return 1;
            default: return 0;
          }

        default: return 0;
      }

    case(7):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (2.592857142857143 + x * (2.6055555555555556 + x * (1.3430555555555554 + x * (0.3888888888888889 + x * (0.06388888888888888 + (0.005555555555555556 + x / 5040.) * x)))));
            case 1: return 1 + x * (1.45 + x * (-0.19444444444444445 + x * (-1.2347222222222223 + x * (-0.7777777777777778 + x * (-0.21388888888888888 + (-0.027777777777777776 - x / 720.) * x)))));
            case 2: return 1 + x * (0.7833333333333333 + x * (-1.05 + x * (-0.9625 + x_sq * (0.175 + (0.05 + x / 240.) * x))));
            case 3: return 1 + x * (0.25 + x * (-1.3611111111111112 + x * (-0.3402777777777778 + x * (0.3888888888888889 + x * (0.09722222222222222 + (-0.027777777777777776 - x / 144.) * x)))));
            case 4: return 1 + x * (-0.25 + x * (-1.3611111111111112 + x * (0.3402777777777778 + x * (0.3888888888888889 + x * (-0.09722222222222222 + (-0.027777777777777776 + x / 144.) * x)))));
            case 5: return 1 + x * (-0.7833333333333333 + x * (-1.05 + x * (0.9625 + x_sq * (-0.175 + (0.05 - x / 240.) * x))));
            case 6: return 1 + x * (-1.45 + x * (-0.19444444444444445 + x * (1.2347222222222223 + x * (-0.7777777777777778 + x * (0.21388888888888888 + (-0.027777777777777776 + x / 720.) * x)))));
            case 7: return 1 + x * (-2.592857142857143 + x * (2.6055555555555556 + x * (-1.3430555555555554 + x * (0.3888888888888889 + x * (-0.06388888888888888 + (0.005555555555555556 - x / 5040.) * x)))));
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 2.592857142857143 + x * (5.211111111111111 + x * (4.029166666666667 + x * (1.5555555555555556 + x * (0.3194444444444444 + (0.03333333333333333 + x / 720.) * x))));
            case 1: return 1.45 + x * (-0.3888888888888889 + x * (-3.7041666666666666 + x * (-3.111111111111111 + x * (-1.0694444444444444 + (-0.16666666666666666 - (7 * x) / 720.) * x))));
            case 2: return 0.7833333333333333 + x * (-2.1 + x * (-2.8875 + x_sq * (0.875 + (0.3 + (7 * x) / 240.) * x)));
            case 3: return 0.25 + x * (-2.7222222222222223 + x * (-1.0208333333333333 + x * (1.5555555555555556 + x * (0.4861111111111111 + (-0.16666666666666666 - (7 * x) / 144.) * x))));
            case 4: return -0.25 + x * (-2.7222222222222223 + x * (1.0208333333333333 + x * (1.5555555555555556 + x * (-0.4861111111111111 + (-0.16666666666666666 + (7 * x) / 144.) * x))));
            case 5: return -0.7833333333333333 + x * (-2.1 + x * (2.8875 + x_sq * (-0.875 + (0.3 - (7 * x) / 240.) * x)));
            case 6: return -1.45 + x * (-0.3888888888888889 + x * (3.7041666666666666 + x * (-3.111111111111111 + x * (1.0694444444444444 + (-0.16666666666666666 + (7 * x) / 720.) * x))));
            case 7: return -2.592857142857143 + x * (5.211111111111111 + x * (-4.029166666666667 + x * (1.5555555555555556 + x * (-0.3194444444444444 + (0.03333333333333333 - x / 720.) * x))));
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 5.211111111111111 + x * (8.058333333333334 + x * (4.666666666666667 + x * (1.2777777777777777 + (0.16666666666666666 + x / 120.) * x)));
            case 1: return -0.3888888888888889 + x * (-7.408333333333333 + x * (-9.333333333333334 + x * (-4.277777777777778 + (-0.8333333333333334 - (7 * x) / 120.) * x)));
            case 2: return -2.1 + x * (-5.775 + x_sq * (3.5 + (1.5 + (7 * x) / 40.) * x));
            case 3: return -2.7222222222222223 + x * (-2.0416666666666665 + x * (4.666666666666667 + x * (1.9444444444444444 + (-0.8333333333333334 - (7 * x) / 24.) * x)));
            case 4: return -2.7222222222222223 + x * (2.0416666666666665 + x * (4.666666666666667 + x * (-1.9444444444444444 + (-0.8333333333333334 + (7 * x) / 24.) * x)));
            case 5: return -2.1 + x * (5.775 + x_sq * (-3.5 + (1.5 - (7 * x) / 40.) * x));
            case 6: return -0.3888888888888889 + x * (7.408333333333333 + x * (-9.333333333333334 + x * (4.277777777777778 + (-0.8333333333333334 + (7 * x) / 120.) * x)));
            case 7: return 5.211111111111111 + x * (-8.058333333333334 + x * (4.666666666666667 + x * (-1.2777777777777777 + (0.16666666666666666 - x / 120.) * x)));
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 8.058333333333334 + x * (9.333333333333334 + x * (3.8333333333333335 + (0.6666666666666666 + x / 24.) * x));
            case 1: return -7.408333333333333 + x * (-18.666666666666668 + x * (-12.833333333333334 + (-3.3333333333333335 - (7 * x) / 24.) * x));
            case 2: return -5.775 + x_sq * (10.5 + (6 + (7 * x) / 8.) * x);
            case 3: return -2.0416666666666665 + x * (9.333333333333334 + x * (5.833333333333333 + (-3.3333333333333335 - (35 * x) / 24.) * x));
            case 4: return 2.0416666666666665 + x * (9.333333333333334 + x * (-5.833333333333333 + x * (-3.3333333333333335 + (35 * x) / 24.)));
            case 5: return 5.775 + x_sq * (-10.5 + (6 - (7 * x) / 8.) * x);
            case 6: return 7.408333333333333 + x * (-18.666666666666668 + x * (12.833333333333334 + (-3.3333333333333335 + (7 * x) / 24.) * x));
            case 7: return -8.058333333333334 + x * (9.333333333333334 + x * (-3.8333333333333335 + (0.6666666666666666 - x / 24.) * x));
            default: return 0;
          }

        case 4:
          switch(id) {
            case 0: return 9.333333333333334 + x * (7.666666666666667 + (2 + x / 6.) * x);
            case 1: return -18.666666666666668 + x * (-25.666666666666668 + (-10 - (7 * x) / 6.) * x);
            case 2: return x * (21 + x * (18 + (7 * x) / 2.));
            case 3: return 9.333333333333334 + x * (11.666666666666666 + (-10 - (35 * x) / 6.) * x);
            case 4: return 9.333333333333334 + x * (-11.666666666666666 + x * (-10 + (35 * x) / 6.));
            case 5: return x * (-21 + (18 - (7 * x) / 2.) * x);
            case 6: return -18.666666666666668 + x * (25.666666666666668 + x * (-10 + (7 * x) / 6.));
            case 7: return 9.333333333333334 + x * (-7.666666666666667 + (2 - x / 6.) * x);
            default: return 0;
          }

        case 5:
          switch(id) {
            case 0: return 7.666666666666667 + (4 + x / 2.) * x;
            case 1: return -25.666666666666668 + (-20 - (7 * x) / 2.) * x;
            case 2: return 21 + x * (36 + (21 * x) / 2.);
            case 3: return 11.666666666666666 + (-20 - (35 * x) / 2.) * x;
            case 4: return -11.666666666666666 + x * (-20 + (35 * x) / 2.);
            case 5: return -21 + (36 - (21 * x) / 2.) * x;
            case 6: return 25.666666666666668 + x * (-20 + (7 * x) / 2.);
            case 7: return -7.666666666666667 + (4 - x / 2.) * x;
            default: return 0;
          }

        case 6:
          switch(id) {
            case 0: return 4 + x;
            case 1: return -20 - 7 * x;
            case 2: return 36 + 21 * x;
            case 3: return -20 - 35 * x;
            case 4: return -20 + 35 * x;
            case 5: return 36 - 21 * x;
            case 6: return -20 + 7 * x;
            case 7: return 4 - x;
            default: return 0;
          }

        case 7:
          switch(id) {
            case 0: return 1;
            case 1: return -7;
            case 2: return 21;
            case 3: return -35;
            case 4: return 35;
            case 5: return -21;
            case 6: return 7;
            case 7: return -1;
            default: return 0;
          }

        default: return 0;
      }

    case(8):
      switch(derivative) {
        case 0:
          switch(id) {
            case 0: return 1 + x * (2.717857142857143 + x * (2.9296626984126983 + x * (1.66875 + x * (0.5567708333333333 + x * (0.1125 + x * (0.013541666666666667 + (0.0008928571428571428 + x / 40320.) * x))))));
            case 1: return 1 + x * (1.5928571428571427 + x * (0.012698412698412698 + x * (-1.2625 + x * (-0.9541666666666667 + x * (-0.325 + x * (-0.058333333333333334 + (-0.005357142857142857 - x / 5040.) * x))))));
            case 2: return 1 + x * (0.95 + x * (-0.9194444444444444 + x * (-1.1375 + x * (-0.16041666666666668 + x * (0.175 + x * (0.07916666666666666 + (0.0125 + x / 1440.) * x))))));
            case 3: return 1 + x * (0.45 + x * (-1.3111111111111111 + x * (-0.6125 + x * (0.32083333333333336 + x * (0.175 + x * (-0.008333333333333333 + (-0.0125 - x / 720.) * x))))));
            case 4: return 1 + x_sq * (-1.4236111111111112 + x_sq * (0.4739583333333333 + x_sq * (-0.052083333333333336 + x_sq / 576.)));
            case 5: return 1 + x * (-0.45 + x * (-1.3111111111111111 + x * (0.6125 + x * (0.32083333333333336 + x * (-0.175 + x * (-0.008333333333333333 + (0.0125 - x / 720.) * x))))));
            case 6: return 1 + x * (-0.95 + x * (-0.9194444444444444 + x * (1.1375 + x * (-0.16041666666666668 + x * (-0.175 + x * (0.07916666666666666 + (-0.0125 + x / 1440.) * x))))));
            case 7: return 1 + x * (-1.5928571428571427 + x * (0.012698412698412698 + x * (1.2625 + x * (-0.9541666666666667 + x * (0.325 + x * (-0.058333333333333334 + (0.005357142857142857 - x / 5040.) * x))))));
            case 8: return 1 + x * (-2.717857142857143 + x * (2.9296626984126983 + x * (-1.66875 + x * (0.5567708333333333 + x * (-0.1125 + x * (0.013541666666666667 + (-0.0008928571428571428 + x / 40320.) * x))))));
            default: return 0;
          }

        case 1:
          switch(id) {
            case 0: return 2.717857142857143 + x * (5.859325396825397 + x * (5.00625 + x * (2.2270833333333333 + x * (0.5625 + x * (0.08125 + (0.00625 + x / 5040.) * x)))));
            case 1: return 1.5928571428571427 + x * (0.025396825396825397 + x * (-3.7875 + x * (-3.816666666666667 + x * (-1.625 + x * (-0.35 + (-0.0375 - x / 630.) * x)))));
            case 2: return 0.95 + x * (-1.8388888888888888 + x * (-3.4125 + x * (-0.6416666666666667 + x * (0.875 + x * (0.475 + (0.0875 + x / 180.) * x)))));
            case 3: return 0.45 + x * (-2.6222222222222222 + x * (-1.8375 + x * (1.2833333333333334 + x * (0.875 + x * (-0.05 + (-0.0875 - x / 90.) * x)))));
            case 4: return x * (-2.8472222222222223 + x_sq * (1.8958333333333333 + x_sq * (-0.3125 + x_sq / 72.)));
            case 5: return -0.45 + x * (-2.6222222222222222 + x * (1.8375 + x * (1.2833333333333334 + x * (-0.875 + x * (-0.05 + (0.0875 - x / 90.) * x)))));
            case 6: return -0.95 + x * (-1.8388888888888888 + x * (3.4125 + x * (-0.6416666666666667 + x * (-0.875 + x * (0.475 + (-0.0875 + x / 180.) * x)))));
            case 7: return -1.5928571428571427 + x * (0.025396825396825397 + x * (3.7875 + x * (-3.816666666666667 + x * (1.625 + x * (-0.35 + (0.0375 - x / 630.) * x)))));
            case 8: return -2.717857142857143 + x * (5.859325396825397 + x * (-5.00625 + x * (2.2270833333333333 + x * (-0.5625 + x * (0.08125 + (-0.00625 + x / 5040.) * x)))));
            default: return 0;
          }

        case 2:
          switch(id) {
            case 0: return 5.859325396825397 + x * (10.0125 + x * (6.68125 + x * (2.25 + x * (0.40625 + (0.0375 + x / 720.) * x))));
            case 1: return 0.025396825396825397 + x * (-7.575 + x * (-11.45 + x * (-6.5 + x * (-1.75 + (-0.225 - x / 90.) * x))));
            case 2: return -1.8388888888888888 + x * (-6.825 + x * (-1.925 + x * (3.5 + x * (2.375 + (0.525 + (7 * x) / 180.) * x))));
            case 3: return -2.6222222222222222 + x * (-3.675 + x * (3.85 + x * (3.5 + x * (-0.25 + (-0.525 - (7 * x) / 90.) * x))));
            case 4: return -2.8472222222222223 + x_sq * (5.6875 + x_sq * (-1.5625 + (7 * x_sq) / 72.));
            case 5: return -2.6222222222222222 + x * (3.675 + x * (3.85 + x * (-3.5 + x * (-0.25 + (0.525 - (7 * x) / 90.) * x))));
            case 6: return -1.8388888888888888 + x * (6.825 + x * (-1.925 + x * (-3.5 + x * (2.375 + (-0.525 + (7 * x) / 180.) * x))));
            case 7: return 0.025396825396825397 + x * (7.575 + x * (-11.45 + x * (6.5 + x * (-1.75 + (0.225 - x / 90.) * x))));
            case 8: return 5.859325396825397 + x * (-10.0125 + x * (6.68125 + x * (-2.25 + x * (0.40625 + (-0.0375 + x / 720.) * x))));
            default: return 0;
          }

        case 3:
          switch(id) {
            case 0: return 10.0125 + x * (13.3625 + x * (6.75 + x * (1.625 + (0.1875 + x / 120.) * x)));
            case 1: return -7.575 + x * (-22.9 + x * (-19.5 + x * (-7 + (-1.125 - x / 15.) * x)));
            case 2: return -6.825 + x * (-3.85 + x * (10.5 + x * (9.5 + (2.625 + (7 * x) / 30.) * x)));
            case 3: return -3.675 + x * (7.7 + x * (10.5 + x * (-1 + (-2.625 - (7 * x) / 15.) * x)));
            case 4: return x * (11.375 + x_sq * (-6.25 + (7 * x_sq) / 12.));
            case 5: return 3.675 + x * (7.7 + x * (-10.5 + x * (-1 + (2.625 - (7 * x) / 15.) * x)));
            case 6: return 6.825 + x * (-3.85 + x * (-10.5 + x * (9.5 + (-2.625 + (7 * x) / 30.) * x)));
            case 7: return 7.575 + x * (-22.9 + x * (19.5 + x * (-7 + (1.125 - x / 15.) * x)));
            case 8: return -10.0125 + x * (13.3625 + x * (-6.75 + x * (1.625 + (-0.1875 + x / 120.) * x)));
            default: return 0;
          }

        case 4:
          switch(id) {
            case 0: return 13.3625 + x * (13.5 + x * (4.875 + (0.75 + x / 24.) * x));
            case 1: return -22.9 + x * (-39 + x * (-21 + (-4.5 - x / 3.) * x));
            case 2: return -3.85 + x * (21 + x * (28.5 + x * (10.5 + (7 * x) / 6.)));
            case 3: return 7.7 + x * (21 + x * (-3 + (-10.5 - (7 * x) / 3.) * x));
            case 4: return 11.375 + x_sq * (-18.75 + (35 * x_sq) / 12.);
            case 5: return 7.7 + x * (-21 + x * (-3 + (10.5 - (7 * x) / 3.) * x));
            case 6: return -3.85 + x * (-21 + x * (28.5 + x * (-10.5 + (7 * x) / 6.)));
            case 7: return -22.9 + x * (39 + x * (-21 + (4.5 - x / 3.) * x));
            case 8: return 13.3625 + x * (-13.5 + x * (4.875 + (-0.75 + x / 24.) * x));
            default: return 0;
          }

        case 5:
          switch(id) {
            case 0: return 13.5 + x * (9.75 + (2.25 + x / 6.) * x);
            case 1: return -39 + x * (-42 + (-13.5 - (4 * x) / 3.) * x);
            case 2: return 21 + x * (57 + x * (31.5 + (14 * x) / 3.));
            case 3: return 21 + x * (-6 + (-31.5 - (28 * x) / 3.) * x);
            case 4: return x * (-37.5 + (35 * x_sq) / 3.);
            case 5: return -21 + x * (-6 + (31.5 - (28 * x) / 3.) * x);
            case 6: return -21 + x * (57 + x * (-31.5 + (14 * x) / 3.));
            case 7: return 39 + x * (-42 + (13.5 - (4 * x) / 3.) * x);
            case 8: return -13.5 + x * (9.75 + (-2.25 + x / 6.) * x);
            default: return 0;
          }

        case 6:
          switch(id) {
            case 0: return 9.75 + (4.5 + x / 2.) * x;
            case 1: return -42 + (-27 - 4 * x) * x;
            case 2: return 57 + x * (63 + 14 * x);
            case 3: return -6 + (-63 - 28 * x) * x;
            case 4: return -37.5 + 35 * x_sq;
            case 5: return -6 + (63 - 28 * x) * x;
            case 6: return 57 + x * (-63 + 14 * x);
            case 7: return -42 + (27 - 4 * x) * x;
            case 8: return 9.75 + (-4.5 + x / 2.) * x;
            default: return 0;
          }

        case 7:
          switch(id) {
            case 0: return 4.5 + x;
            case 1: return -27 - 8 * x;
            case 2: return 63 + 28 * x;
            case 3: return -63 - 56 * x;
            case 4: return 70 * x;
            case 5: return 63 - 56 * x;
            case 6: return -63 + 28 * x;
            case 7: return 27 - 8 * x;
            case 8: return -4.5 + x;
            default: return 0;
          }

        case 8:
          switch(id) {
            case 0: return 1;
            case 1: return -8;
            case 2: return 28;
            case 3: return -56;
            case 4: return 70;
            case 5: return -56;
            case 6: return 28;
            case 7: return -8;
            case 8: return 1;
            default: return 0;
          }

        default: return 0;
      }

    default: return 0;
  }
}
//clang-format on
