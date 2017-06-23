#ifndef CMPLX_ODE_RHS
#define CMPLX_ODE_RHS


namespace Integrator {
  class CMPLX_0DE_RHS;
}

class Integrator::CMPLX_ODE_RHS : public RHS<std::complex<double>> {
  public:
   CMPLX_ODE_RHS(const double dt, const std::shared_ptr<
		 History<std::complex<double>> &history)
       : dt(dt), history(history){};
   void evaluate(const int) const;

  private:
};

void Integrator::CMPLX_ODE_RHS::evaluate(const int n) const
{
  const std::complex<double> iu(0,1);
  const double time = n * dt;
  for(int i; i<static_cast<int>(history->array.shape()[0]); ++i) {
    history->array[i][n][1] = 1;
  }
}

#endif
	
