#ifndef JASTROW_POSCH_TELLER_H
#define JASTROW_POSCH_TELLER_H


#include "jastrow.h"

class jastrowPoschTeller  : public jastrow<jastrowPoschTeller>
{
public:
  
  jastrowPoschTeller(real_t R0_,real_t C_,real_t alpha_,real_t Rm_,real_t lBox_);
  jastrowPoschTeller(const json_t & j);
  
  static std::string name() {return "poschTeller";}
  
  void initCoefficients();
  
  inline double d0(const double & r) const {return r<=Rm ? d0_1(r) :
                                             (   
                                            
					     r<=halflBox ? d0_2(r) : longDistanceConstant
						);
					  
                                             }

  inline double d1(const double & r) const {return
                                              r<= Rm ? d1_1(r) :
					      (
					       r<=halflBox ? d1_2(r) : 0
						  );
						 
                                           }
  
  inline double d2(const double & r) const {return   
                                              r<= Rm ? d2_1(r) :
					      (
					         r<=halflBox ? d2_2(r) : 0
					      );
						 
                                           }
  
  inline double d0_1(const double & r) const {return log( tanh(K0*r)/r);}
  inline double d1_1(const double & r) const {return (2*K0/sinh(2*K0*r) -1/(r));}
  inline double d2_1(const double & r) const {return (-4*K0*K0*cosh(2*K0*r)/(sinh(2*K0*r)*sinh(2*K0*r) ) +1./(r*r));}
  
  
  inline double d0_2(const double & r) const {return log(C*(exp(-alpha*r) + exp(-alpha*(lBox-r))));};
  inline double d1_2(const double & r) const {return alpha *(-exp(-alpha*r) + exp(-alpha*(lBox-r))) /   ((exp(-alpha*r) + exp(-alpha*(lBox-r))));   }
  inline double d2_2(const double & r) const {
    double tmp=(exp(-alpha*r) + exp(-alpha*(lBox-r)));
    double tmp1=(-exp(-alpha*r) + exp(-alpha*(lBox-r)));
    return alpha*alpha*( 1. -   tmp1*tmp1/pow(tmp,2)    );
  }
  
private:
  
  double alpha;
  double K0;
  double R0;  
  double C;
  double longDistanceConstant;
  double lBox;
  double Rm;
  double halflBox;
};


#endif
