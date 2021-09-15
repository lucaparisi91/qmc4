#include "jastrow.h"
#include <cmath>


class jastrowSquareWell  : public jastrow<jastrowSquareWell>
{
public:
  
  jastrowSquareWell(real_t V0_,real_t R0_,real_t Rm_,real_t alpha_,real_t aInverse_,real_t lBox_) : V0(V0_),R0(R0_),Rm(Rm_),alpha(alpha_),lBox(lBox_),aInverse(aInverse_)
  {
    initCoefficients();
  }
  
  void initCoefficients();

  jastrowSquareWell(const json_t & j) ;
  
  static std::string name() {return "squareWell";}
  
  inline real_t d0(const real_t & r) const {return r<=R0 ? logd0_1(r) :
                                             (   
                                              r<= Rm ? logd0_2(r) :
					      (
					         r<=halflBox ? logd0_3(r) : longDistanceConstant
					       )
						 );
                                           }

  inline real_t d1(const real_t & r) const {return r<=R0 ? logd1_1(r) :
      (   
                                              r<= Rm ? logd1_2(r) :
					      (
					       r<=halflBox ? logd1_3(r) : 0
					       )
						 );
                                           }
  
  inline real_t d2(const real_t & r) const {return r<=R0 ? logd2_1(r) :
                                             (   
                                              r<= Rm ? logd2_2(r) :
					      (
					         r<=halflBox ? logd2_3(r) : 0
					       )
						 );
                                           }
private:
  inline real_t logd0_1(const real_t & r) const {return log( std::sin(K0*r)/r);}
  inline real_t logd1_1(const real_t & r) const {return (K0/std::tan(K0*r)-1/r);}
  inline real_t logd2_1(const real_t & r) const {return (-K0*K0/std::pow(sin(K0*r),2) +1/(r*r));}
  
  inline real_t logd0_2(const real_t & r) const {return std::log( 1/r - aInverse );};
  inline real_t logd1_2(const real_t & r) const {return 1/(r*(-1 + r*aInverse));}
  inline real_t logd2_2(const real_t & r) const  {return (-2*r*aInverse+1)/(r*r*(r*aInverse-1)*(r*aInverse-1));}
  
  
  inline real_t logd0_3(const real_t & r) const {return std::log(B+C*(std::exp(-alpha*r) + std::exp(-alpha*(lBox-r))));};
  inline real_t logd1_3(const real_t & r) const {return alpha * C*(-std::exp(-alpha*r) + std::exp(-alpha*(lBox-r))) /   (B+C*(std::exp(-alpha*r) + std::exp(-alpha*(lBox-r))));   }
  inline real_t logd2_3(const real_t & r) const {
  real_t tmp=(exp(-alpha*r) + exp(-alpha*(lBox-r)));
  real_t tmp1=(-exp(-alpha*r) + exp(-alpha*(lBox-r)));
  return alpha*alpha*C *( (tmp)/(B+C*tmp) -   C*tmp1*tmp1/pow(B+C*tmp,2)    );
  }
  
  real_t aInverse;//inverse scattering length
  real_t R0;//range of the potential
  real_t Rm; // matching point for the potential
  real_t alpha;
  real_t V0;
  
  real_t B,C,A;
  real_t K0;
  real_t longDistanceConstant;
  real_t lBox;
  real_t halflBox;
};


class jastrowSquareWellBoundState  : public jastrow<jastrowSquareWellBoundState>
{
public:
  static std::string name() {return "squareWellBoundState";}

  jastrowSquareWellBoundState(real_t V0_,real_t R0_,real_t Rm_,real_t lBox_,real_t _aInverse,real_t _eb,real_t alpha_) : V0(V0_),R0(R0_),Rm(Rm_),lBox(lBox_),aInverse(_aInverse),eb(_eb),alpha(alpha_)
  {
    initCoefficients();
  }

  jastrowSquareWellBoundState(const json_t & j);
  
  void initCoefficients();
  
  
  inline real_t d0(const real_t & r) const {return r<=R0 ? d0_1(r) :
                                             (   
                                              r<= Rm ? d0_2(r) :
					      (
					         r<=halflBox ? d0_3(r) : longDistanceConstant
					       )
						 );
                                           }

  inline real_t d1(const real_t & r) const {return r<=R0 ? d1_1(r) :
      (   
                                              r<= Rm ? d1_2(r) :
					      (
					       r<=halflBox ? d1_3(r) : 0
					       )
						 );
                                           }
  
  inline real_t d2(const real_t & r) const {return r<=R0 ? d2_1(r) :
                                             (   
                                              r<= Rm ? d2_2(r) :
					      (
					         r<=halflBox ? d2_3(r) : 0
					       )
						 );
                                           }
  
  inline real_t d0_1(const real_t & r) const {return log( sin(K0*r)/r) + logA;}
  inline real_t d1_1(const real_t & r) const {return (K0/tan(K0*r)-1/r);}
  inline real_t d2_1(const real_t & r) const {return (-K0*K0/pow(sin(K0*r),2) +1/(r*r));}
  
  inline real_t d0_2(const real_t & r) const {return - expCoeff*r - log(r) ;}
  inline real_t d1_2(const real_t & r) const {return -expCoeff -1./r;}
  inline real_t d2_2(const real_t & r) const {return 1./(r*r);}
  
  inline real_t d0_3(const real_t & r) const {return log( B + C*(exp(-alpha*r) + exp(-alpha*(lBox-r))));};
  inline real_t d1_3(const real_t & r) const {return alpha * C*(-exp(-alpha*r) + exp(-alpha*(lBox-r))) /   (B + C*(exp(-alpha*r) + exp(-alpha*(lBox-r))));   }
  inline real_t d2_3(const real_t & r) const{
    real_t tmp=(exp(-alpha*r) + exp(-alpha*(lBox-r)));
    real_t tmp1=(-exp(-alpha*r) + exp(-alpha*(lBox-r)));
    return alpha*alpha*C *( (tmp)/(B+ C*tmp) -   C*tmp1*tmp1/pow(B + C*tmp,2)   );
  }
  
  
private:
  
  real_t aInverse;//inverse scattering length
  real_t R0;//range of the potential
  real_t Rm; // matching point for the potential
  real_t alpha;
  real_t V0,eb;
  real_t m;
  real_t C,A,B,logA,expCoeff;
  real_t K0;
  real_t longDistanceConstant;
  real_t lBox;
  real_t halflBox;
  
};
