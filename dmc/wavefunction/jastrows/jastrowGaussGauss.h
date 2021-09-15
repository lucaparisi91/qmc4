#ifndef JASTROWGAUSSGAUSS_H
#define JASTROWGAUSSGAUSS_H


#include "jastrow.h"


class jastrowGaussGauss: public jastrow<jastrowGaussGauss>
{
	// J(x) = -alpha* x^2 
public:
  
  jastrowGaussGauss(real_t r0_,real_t u0_);
  jastrowGaussGauss(const json_t  & j);

  
  static std::string name() {return "gaussiangaussian"; }
  
  real_t d0(real_t x) const {return u0*exp(-x*x*alpha);}
  real_t d1(real_t x) const {return -2*u0*alpha*x*exp(-x*x*alpha);}
  real_t d2(real_t x) const {return -2*alpha*u0*(1 - 2*alpha*x*x)*exp(-x*x*alpha);}
  
private:
  real_t alpha;
  real_t u0;
  
  
};


#endif
