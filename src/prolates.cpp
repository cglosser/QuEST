#include "prolates.h"

namespace math = boost::math::constants;

double sinc(double);
double sinhc(double);

double d1_sinc(double);
double d1_sinhc(double);

double d2_sinc(double);
double d2_sinhc(double);

double d3_sinc(double);
double d3_sinhc(double);

Prolate::Prolate(int n)
{
  width = n;
  alpha = width*math::pi<double>();
}

double Prolate::sqrt_term(const double t)  {return std::sqrt(1 - std::pow(t/width, 2));}
double Prolate::d1_sqrt_term(const double t) {return -t/(sqrt_term(t)*pow(width, 2));}
double Prolate::d2_sqrt_term(const double t)
{
  const double q = sqrt_term(t);
  return -std::pow(t,2)/(width*std::pow(q*width,3)) - 1/(q*std::pow(width,2));
}
double Prolate::d3_sqrt_term(const double t)
{
  const double q = sqrt_term(t);
  return -3*t/width*(std::pow(t,2)/std::pow(q*width,5) + 1/std::pow(q*width,3));
}

//=============================================================

double Prolate::d0(const double t)
{
  if(std::abs(t) >= width) return 0;
  const double pi = math::pi<double>(), q = sqrt_term(t);
  return sinc(pi*t)*alpha*sinhc(alpha*q)/std::sinh(alpha);
}

double Prolate::d1(const double t)
{
  const double pi = math::pi<double>();

  if(std::abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t);
    return
      ( //If you look closely, you'll see a product rule here.
        alpha*pi*sinhc(alpha*q)*d1_sinc(pi*t) +
        std::pow(alpha, 2)*d1q*sinc(pi*t)*d1_sinhc(alpha*q)
      )/std::sinh(alpha);
  }
}

double Prolate::d2(const double t)
{
  const double pi = math::pi<double>();

  if(std::abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t), d2q = d2_sqrt_term(t);
    return
      (
        2*std::pow(alpha,2)*d1q*pi*d1_sinc(pi*t)*d1_sinhc(alpha*q) +
        alpha*std::pow(pi,2)*sinhc(alpha*q)*d2_sinc(pi*t) +
        alpha*sinc(pi*t)*(alpha*d2q*d1_sinhc(alpha*q) + std::pow(alpha*d1q,2)*d2_sinhc(alpha*q))
      )/std::sinh(alpha);
  }
}

double Prolate::d3(const double t)
{
  const double pi = math::pi<double>();

  if (std::abs(t) >= width) {
    return 0;
  } else {
    const double   q =    sqrt_term(t),  d2q = d2_sqrt_term(t);
    const double d1q = d1_sqrt_term(t) , d3q = d3_sqrt_term(t);

    double x;
    std::cout << q << std::endl;
    std::cout << d1q << std::endl;
    std::cout << d2q << std::endl;
    std::cout << d3q << std::endl;

    std::cin >> x;

    return (
alpha*(3*alpha*d1q*std::pow(pi,2)*d1_sinhc(alpha*q)*d2_sinc(pi*t) +
     3*alpha*pi*d1_sinc(pi*t)*(d2q*d1_sinhc(alpha*q) + alpha*std::pow(d1q,2)*d2_sinhc(alpha*q)) +
     alpha*(d3q*d1_sinhc(alpha*q) + alpha*d1q*
         (3*d2q*d2_sinhc(alpha*q) + alpha*std::pow(d1q,2)*d3_sinhc(alpha*q)))*sinc(pi*t) +
     std::pow(pi,3)*d3_sinc(pi*t)*sinhc(alpha*q))
      )/std::sinh(alpha);
  }
}

//=============================================================

double sinc(const double t)
{
  if (abs(t) <= TOLER) {
    return 1 - std::pow(t, 2)/6 + std::pow(t, 4)/120;
  } else {
    return sin(t)/t;
  }
}

double sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return 1 + std::pow(t, 2)/6 + std::pow(t, 4)/120;
  } else {
    return std::sinh(t)/t;
  }
}

//=============================================================

double d1_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return -t/3 + std::pow(t, 3)/30 - std::pow(t, 5)/840;
  } else {
    return cos(t)/t - sin(t)/std::pow(t, 2);
  }
}

double d1_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/3 + std::pow(t, 3)/30 + std::pow(t, 5)/840;
  } else {
    return cosh(t)/t - sinh(t)/std::pow(t, 2);
  }
}

//=============================================================

double d2_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return -1.0/3 + std::pow(t,2)/10 - std::pow(t,4)/168;
  } else {
    return -std::sin(t)/t - 2*std::cos(t)/std::pow(t,2) + 2*std::sin(t)/pow(t,3);
  }
}

double d2_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return 1.0/3 + std::pow(t,2)/10 + std::pow(t,4)/168;
  } else {
    return std::sinh(t)/t - 2*std::cosh(t)/std::pow(t,2) + 2*std::sinh(t)/pow(t,3);
  }
}

//=============================================================

double d3_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/5 - std::pow(t,3)/42 + std::pow(t,5)/1080;
  } else {
    return -std::cos(t)/t + 3*std::sin(t)/std::pow(t,2) +
            6*std::cos(t)/std::pow(t,3) - 6*std::sin(t)/std::pow(t,4);
  }
}

double d3_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/5 + std::pow(t,3)/42 + std::pow(t,5)/1080;
  } else {
    return std::cos(t)/t - 3*std::sin(t)/std::pow(t,2) +
            6*std::cos(t)/std::pow(t,3) - 6*std::sin(t)/std::pow(t,4);
  }
}
