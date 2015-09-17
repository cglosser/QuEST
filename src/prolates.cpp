#include "prolates.h"

using namespace std;

Prolate::Prolate(double w0, double timestep)
{
  omega0 = w0;
  dt = timestep;
}

//=============================================================

double sinc(double t, double omega = 1)
{
  double arg = omega*t;
  if (abs(arg) <= TOLER) {
    return 1 - pow(arg, 2)/6 + pow(arg, 4)/120;
  } else {
    return sin(arg)/arg;
  }
}

double dsinc(double t, double omega = 1)
{
  double arg = omega*t;
  if (abs(arg) <= TOLER) {
    return omega*(0 - arg/3 + pow(arg, 3)/30);
  } else {
    return omega*(cos(arg)/arg - sin(arg)/pow(arg, 2));
  }
}

double ddsinc(double t, double omega = 1)
{
  double arg = omega*t;
  if (abs(arg) <= TOLER) {
    return -1*pow(omega, 2)/3*(1 - 3*pow(arg, 2)/10 + pow(arg, 4)/56);
  } else {
    return pow(omega, 2)*
      (-2*cos(arg)/pow(arg, 2) - sin(arg)/arg + 2*sin(arg)/pow(arg, 3));
  }
}

double dddsinc(double t, double omega = 1)
{
  double arg = omega*t;
  if (abs(arg) <= TOLER) {
    return 0 + t*pow(omega, 4)/5 - pow(t, 3)*pow(omega, 6)/42;
  } else {
    return pow(omega, 3)*(
      6*cos(arg)/pow(arg, 3) -   cos(arg)/arg - 
      6*sin(arg)/pow(arg, 4) + 3*sin(arg)/pow(arg, 2)
    );
  }
}

//=============================================================

double sinhc(double t, double omega = 1)
{
  double arg = omega*t;
  if (abs(arg) <= TOLER) {
    return 1 + pow(arg, 2)/6 + pow(arg, 4)/120;
  } else {
    return sinh(arg)/arg;
  }
}
