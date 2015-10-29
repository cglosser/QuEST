#include "prolates.h"

constexpr double pi = boost::math::constants::pi<double>();

using std::abs;
using std::cos;
using std::cosh;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;

double sinc(const double);
double sinhc(const double);

double d1_sinc(const double);
double d1_sinhc(const double);

double d2_sinc(const double);
double d2_sinhc(const double);

double d3_sinc(const double);
double d3_sinhc(const double);

Prolate::Prolate(const int n)
{
  width = n;
  alpha = width*pi;
}

double Prolate::sqrt_term(const double t) const
{
  return sqrt(1 - pow(t/width, 2));
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
  if(abs(t) >= width) return 0;
  const double q = sqrt_term(t);
  return sinc(pi*t)*alpha*sinhc(alpha*q)/sinh(alpha);
}

double Prolate::d1(const double t) const
{
  if(abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t);
    return (
        //If you look closely, you'll see a product rule here.
        alpha*pi*sinhc(alpha*q)*d1_sinc(pi*t) +
        pow(alpha, 2)*d1q*sinc(pi*t)*d1_sinhc(alpha*q)
      )/sinh(alpha);
  }
}

double Prolate::d2(const double t) const
{
  if(abs(t) >= width) {
    return 0;
  } else {
    const double q = sqrt_term(t), d1q = d1_sqrt_term(t), d2q = d2_sqrt_term(t);
    return (
        2*pow(alpha,2)*d1q*pi*d1_sinc(pi*t)*d1_sinhc(alpha*q) +
        alpha*pow(pi,2)*sinhc(alpha*q)*d2_sinc(pi*t) +
        alpha*sinc(pi*t)*(alpha*d2q*d1_sinhc(alpha*q) + pow(alpha*d1q,2)*d2_sinhc(alpha*q))
      )/sinh(alpha);
  }
}

{
  }
}

//=============================================================

double sinc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return 1 + t_sq*(-1.0/6 + t_sq/120); 
  } else {
    return sin(t)/t;
  }
}

double sinhc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return 1 + t_sq*(1.0/6 + t_sq/120);
  } else {
    return sinh(t)/t;
  }
}

//=============================================================

double d1_sinc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return t*(-1.0/3 + t_sq*(1.0/30 - t_sq/840));
  } else {
    return cos(t)/t - sin(t)/pow(t, 2);
  }
}

double d1_sinhc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return t*(1.0/3 + t_sq*(1.0/30 + t_sq/840));
  } else {
    return cosh(t)/t - sinh(t)/pow(t, 2);
  }
}

//=============================================================

double d2_sinc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return -1.0/3 + t_sq*(1.0/10 - t_sq/168);
  } else {
    return -sin(t)/t - 2*cos(t)/pow(t,2) + 2*sin(t)/pow(t,3);
  }
}

double d2_sinhc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return 1.0/3 + t_sq*(1.0/10 + t_sq/168);
  } else {
    return sinh(t)/t - 2*cosh(t)/pow(t,2) + 2*sinh(t)/pow(t,3);
  }
}

//=============================================================

double d3_sinc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return t*(1.0/5 + t_sq*(-1.0/42 + t_sq/1080));
  } else {
    return -cos(t)/t + 3*sin(t)/pow(t,2) +
            6*cos(t)/pow(t,3) - 6*sin(t)/pow(t,4);
  }
}

double d3_sinhc(const double t)
{
  if (abs(t) <= TOLER) {
    double t_sq = pow(t, 2);
    return t*(1.0/5 + t_sq*(1.0/42 + t_sq/1080));
  } else {
    return cos(t)/t - 3*sin(t)/pow(t,2) +
            6*cos(t)/pow(t,3) - 6*sin(t)/pow(t,4);
  }
}
