#include "jastrow_delta.h"
#include <nlohmann/json.hpp>


jastrow_delta_phonons::jastrow_delta_phonons(const json_t & j)
{

  parameters.resize(6);
  
  parameters[0]=j["k"].get<real_t>();
  parameters[1]=j["delta"].get<real_t>();
  parameters[2]=j["beta"].get<real_t>();
  parameters[3]=j["z"].get<real_t>();
  parameters[4]=j["cut_off"].get<real_t>();
  k2=M_PI/(2*parameters[4]);

  try
  {
    parameters[5]=j["g"].get<real_t>();
  
  }
  catch(invalidType & e)
  {    
    if (j["g"].get<std::string>() == "inf"  )
    {
      parameters[5]=0;
    }
    else
    {
      throw invalidInput("g is not a number and is not infinite");
    }

    
  }
  

  
  
}


void jastrow_delta_phonons::process()
 {
    double aP=1e-5;
    double bP=3*aP;
    
    fR.z=parameters[3];
    fR.g=parameters[5];
    fR.l_box=parameters[4]*2;
    
    while( fR(bP)*fR(aP) >=0 )
    {

      
      bP+=fR.l_box*1e-4;
      if (bP > fR.l_box)
       	{
	  
	  std::cout << "Failed to find a root. " << std::endl;
	  exit(1);
	}
    }
        parameters[0]=findRootBrente(fR,aP,bP,1e-12);
    parameters[1]=atan(parameters[0]/parameters[5]);

    parameters[2]=(parameters[0]*tan((M_PI/(2*parameters[4]))*parameters[3]))/(M_PI/(parameters[4]*2)  * tan(parameters[0]*parameters[3]+parameters[1]));
    
    k2=M_PI/(2*parameters[4]);
    
    
  }


jastrow_delta_in_trap::jastrow_delta_in_trap(const json_t & j) : jastrow_delta_in_trap::jastrow_delta_in_trap(j["a"].get<real_t>())
{
  
}

jastrow_delta_bound_state_phonons::jastrow_delta_bound_state_phonons(const json_t & j)
{  
  xI=j["xI"].get<real_t>();
  g=j["g"].get<real_t>();  
  lBox=j["cut_off"].get<real_t>()*2;
  k=j["k"].get<real_t>();
  A=j["A"].get<real_t>();
  logA=log(A);
  beta=j["beta"].get<real_t>();
  k2=M_PI/lBox;
}
