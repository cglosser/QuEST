#include "prolates.h"

using boost::math::constants::pi;
using std::pow;

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
  alpha = width*pi<double>();
}

double Prolate::sqrt_term(const double t) const
{
  return std::sqrt(1 - pow(t/width, 2));
}

double Prolate::d1_sqrt_term(const double t) const
{
  return -t/(sqrt_term(t)*pow(width, 2));
}

double Prolate::d2_sqrt_term(const double t) const
{
  const double q = sqrt_term(t);
  return -pow(t,2)/(width*pow(q*width,3)) - 1/(q*pow(width,2));
}

double Prolate::d3_sqrt_term(const double t) const
{
  const double q = sqrt_term(t);
  return -3*t/width*(pow(t,2)/pow(q*width,5) + 1/pow(q*width,3));
}

//=============================================================

double Prolate::d0(const double t) const
{
  if(std::abs(t) >= width) return 0;
  const double q = sqrt_term(t);
  return sinc(pi<double>()*t)*alpha*sinhc(alpha*q)/std::sinh(alpha);
}

double Prolate::d1(const double t) const
{
  if(std::abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t);
    return (
        //If you look closely, you'll see a product rule here.
        alpha*pi<double>()*sinhc(alpha*q)*d1_sinc(pi<double>()*t) +
        pow(alpha, 2)*d1q*sinc(pi<double>()*t)*d1_sinhc(alpha*q)
      )/std::sinh(alpha);
  }
}

double Prolate::d2(const double t) const
{
  if(std::abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t), d2q = d2_sqrt_term(t);
    return (
        2*pow(alpha,2)*d1q*pi<double>()*d1_sinc(pi<double>()*t)*d1_sinhc(alpha*q) +
        alpha*pow(pi<double>(),2)*sinhc(alpha*q)*d2_sinc(pi<double>()*t) +
        alpha*sinc(pi<double>()*t)*(alpha*d2q*d1_sinhc(alpha*q) + pow(alpha*d1q,2)*d2_sinhc(alpha*q))
      )/std::sinh(alpha);
  }
}

double Prolate::d3(const double t) const
{
  if (std::abs(t) >= width) {
    return 0;
  } else {
    const double   q =    sqrt_term(t),  d2q = d2_sqrt_term(t);
    const double d1q = d1_sqrt_term(t) , d3q = d3_sqrt_term(t);
    return (
        3*alpha*d1q*pow(pi<double>(),2)*d1_sinhc(alpha*q)*d2_sinc(pi<double>()*t) +
        3*alpha*pi<double>()*d1_sinc(pi<double>()*t)*(d2q*d1_sinhc(alpha*q) + alpha*pow(d1q,2)*d2_sinhc(alpha*q)) +
        alpha*(d3q*d1_sinhc(alpha*q) + alpha*d1q*(3*d2q*d2_sinhc(alpha*q) + alpha*pow(d1q,2)*d3_sinhc(alpha*q)))*sinc(pi<double>()*t) +
        pow(pi<double>(),3)*d3_sinc(pi<double>()*t)*sinhc(alpha*q)
      )*alpha/std::sinh(alpha);
  }
}

//=============================================================

double sinc(const double t)
{
  if (abs(t) <= TOLER) {
    return 1 - pow(t, 2)/6 + pow(t, 4)/120;
  } else {
    return sin(t)/t;
  }
}

double sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return 1 + pow(t, 2)/6 + pow(t, 4)/120;
  } else {
    return std::sinh(t)/t;
  }
}

//=============================================================

double d1_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return -t/3 + pow(t, 3)/30 - pow(t, 5)/840;
  } else {
    return cos(t)/t - sin(t)/pow(t, 2);
  }
}

double d1_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/3 + pow(t, 3)/30 + pow(t, 5)/840;
  } else {
    return cosh(t)/t - sinh(t)/pow(t, 2);
  }
}

//=============================================================

double d2_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return -1.0/3 + pow(t,2)/10 - pow(t,4)/168;
  } else {
    return -std::sin(t)/t - 2*std::cos(t)/pow(t,2) + 2*std::sin(t)/pow(t,3);
  }
}

double d2_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return 1.0/3 + pow(t,2)/10 + pow(t,4)/168;
  } else {
    return std::sinh(t)/t - 2*std::cosh(t)/pow(t,2) + 2*std::sinh(t)/pow(t,3);
  }
}

//=============================================================

double d3_sinc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/5 - pow(t,3)/42 + pow(t,5)/1080;
  } else {
    return -std::cos(t)/t + 3*std::sin(t)/pow(t,2) +
            6*std::cos(t)/pow(t,3) - 6*std::sin(t)/pow(t,4);
  }
}

double d3_sinhc(const double t)
{
  if (std::abs(t) <= TOLER) {
    return t/5 + pow(t,3)/42 + pow(t,5)/1080;
  } else {
    return std::cos(t)/t - 3*std::sin(t)/pow(t,2) +
            6*std::cos(t)/pow(t,3) - 6*std::sin(t)/pow(t,4);
  }
}
