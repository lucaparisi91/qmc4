#include "jastrowGaussGauss.h"



jastrowGaussGauss::jastrowGaussGauss(real_t r0,real_t u0_) :
  alpha(1./(2*r0*r0)),
  u0(u0_)
{
  
}

jastrowGaussGauss::jastrowGaussGauss(const json_t & j) :
  jastrowGaussGauss::jastrowGaussGauss(j["r0"].get<real_t>(),j["u0"].get<real_t>() )
{
  
  
};


