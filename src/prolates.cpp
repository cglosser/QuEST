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

//=============================================================
//=====  Prolate basis functions                          =====
//=============================================================

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

//=============================================================
//=====  Prolate Time Expansions                          =====
//=============================================================

ProlateTimeExpansion::ProlateTimeExpansion(const Prolate &bf, const int duration,
    const Eigen::Vector3d &init)
{
  basis_function = bf;
  history.reserve(duration);
  history.push_back(init);
}

//=============================================================

void ProlateTimeExpansion::step(const Eigen::Vector3d &new_val)
{
  history.push_back(new_val);
}

//=============================================================

Eigen::Vector3d ProlateTimeExpansion::at(const double time)
{
  /* To help modularity, the time expansion works in terms of *dimensionless*
   * times (t/dt). Thus you must also make the argument of this function
   * dimensionless for an appropriate query of the signal's value. */

  //assert(time <= history.size()); //Can't look into the future

  double intpart, fracpart;
  fracpart = std::modf(time, &intpart);

  int lbound = std::floor(intpart) - basis_function.get_width() + 1;
  int ubound = std::floor(intpart) + basis_function.get_width();

  Eigen::Vector3d signal(0,0,0);

  if(lbound < 0) {
    for(int i = lbound; i < 0; ++i) {
      // Assume constant initial value for time <= 0
      signal += basis_function.d0(-i + fracpart)*history[0];
    }
    for(int i = 0; i <= ubound; ++i) {
      signal += basis_function.d0(i - time)*history[i];
    }
  } else if ((size_t)ubound > history.size()) {
    assert(false && "Disfuturification not implemented!");
  } else {
    for(int i = lbound; i <= ubound; ++i) {
      signal += basis_function.d0(i - time)*history[i];
    }
  }

  return signal;
}

//=============================================================
//=====  Helper functions                                 =====
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

