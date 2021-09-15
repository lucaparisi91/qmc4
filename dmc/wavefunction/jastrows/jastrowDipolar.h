#ifndef JASTROW_DIPOLAR_H
#define JASTROW_DIPOLAR_H

#include "jastrow.h"
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/bessel_prime.hpp>


class jastrowDipolarRep : public jastrow<jastrowDipolarRep>
{
public:
  
  jastrowDipolarRep(const json_t & j);
  
  jastrowDipolarRep(double D_,double lBox_,double matching_point_,double alpha_,double C_) : D(D_),lBox(lBox_),alpha(alpha_),C(C_),matching_point(matching_point_)
  {
     evaluateParameters();
  }
  
  inline double d0(const double & x) const {return x  < matching_point ? logd0Dipolar(x) + logC : ( x<=lBox/2. ? logd0Tail(x)  : 0) ;}
  inline double d1(const double & x) const {return x  < matching_point ? logd1Dipolar(x)  : ( x<=lBox/2. ? logd1Tail(x)  : 0) ;}
  inline double d2(const double & x) const {return x  < matching_point ? logd2Dipolar(x) : ( x<=lBox/2. ? logd2Tail(x)  : 0) ;}
  
  static std::string name() {return "dipolar_rep";}
  
private:
  
  void evaluateParameters();

  
  inline double logd0Dipolar(double x) const {double a = sqrt(x);return log(a* boost::math::cyl_bessel_k(1, 2*D/a)) ; }
  
  inline double logd1Dipolar(double x) const {
    double a = sqrt(x);
    double y=2*D/a;
    double J=boost::math::cyl_bessel_k(1, y);
    double tmp1=boost::math::cyl_bessel_k(0, y);
    double tmp2=boost::math::cyl_bessel_k(2, y);
    double Jprime=-(tmp2 + tmp1)/2.;
    Jprime/=J;
    
    return 0.5/x*(1 -   Jprime*y) ; }
  
  inline double logd2Dipolar(double x) const {
    // compute bessel functions and their inverses
    double a=sqrt(x);
    double y=2*D/a;
    double J=boost::math::cyl_bessel_k(1, y);
    double tmp1=boost::math::cyl_bessel_k(0, y);
    double tmp2=boost::math::cyl_bessel_k(2, y);
    double tmp3=boost::math::cyl_bessel_k(3, y);    
    double tmp4=boost::math::cyl_bessel_k(-1, y);
    
    double Jprime=-(tmp2 + tmp1)/2.;
    double Jprime2=(tmp3 + tmp4 + 2*J)/4.;
    
    Jprime=Jprime/J;
    Jprime2=Jprime2/J;
    
    return ( -1 + y*Jprime*3/2. + 0.5*y*y*(Jprime2-Jprime*Jprime))/(2*x*x);
    
  }
  
  inline double logd0Tail(double x) const  {return alpha*log(sin(k*x));}
  inline double logd1Tail(double x) const { return alpha*k*1./tan(k*x);}
  inline double logd2Tail(double x) const {return -alpha*k*k/std::pow(sin(k*x),2);};
  
  
  
  double k;
  double C;
  double D;
  double lBox;
  double alpha;
  double matching_point;
  double logC;
};


#endif
