#include "jastrowPoschTeller.h"


void jastrowPoschTeller::initCoefficients()
{
  K0=1/R0;
  halflBox=lBox/2.;
  longDistanceConstant=2*C*exp(-alpha*halflBox);
  
}

jastrowPoschTeller::jastrowPoschTeller(real_t R0_,real_t C_,real_t alpha_,real_t Rm_,real_t lBox_) : R0(R0_),C(C_),alpha(alpha_),Rm(Rm_),lBox(lBox_)
{  
  initCoefficients();
}

jastrowPoschTeller::jastrowPoschTeller(const json_t & j)
  : jastrowPoschTeller (j["R0"].get<real_t>() , j["C"].get<real_t>(),j["alpha"].get<real_t>() ,j["Rm"].get<real_t>() ,j["cut_off"].get<real_t>()*2.)
    
{}
