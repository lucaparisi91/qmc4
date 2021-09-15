#include "jastrowDipolar.h"

jastrowDipolarRep::jastrowDipolarRep(const json_t & j) :
  jastrowDipolarRep(
		    sqrt(j["D"].get<real_t>()),
		    j["cut_off"].get<real_t>()*2,
		    j["matching_point"].get<real_t>(),
		    j["alpha"].get<real_t>(),
		    j["C"].get<real_t>()
		    )
{
  
  
 
  
}

void jastrowDipolarRep::evaluateParameters()
{  
  k=M_PI/lBox;
  
  logC=log(C);
}
