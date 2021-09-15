#include "jastrow.h"
#include "tools.h"

class jastrow_delta_phonons : public jastrow<jastrow_delta_phonons>
{
  typedef double position_t;
  typedef real_t value_t;
public:
  
  static std::string name() {return "delta_phonons";}

  jastrow_delta_phonons(const json_t & j);
  
  
  inline real_t  d0(const double & x) const {return log (d0Exp(x));}
  
  inline real_t d1(const double & x) const {return (x < parameters[3]) ? parameters[0]/tan(parameters[0]*x + parameters[1]) : x > parameters[4] ? 0 : k2*parameters[2]/tan(k2*x) ; }
  inline real_t d2(const double & x) const {return (x < parameters[3]) ? -parameters[0]*parameters[0]/pow(sin(parameters[0]*x + parameters[1]),2) : x > parameters[4] ? 0 : -k2*k2*parameters[2]/pow(sin(k2*x),2) ; }
  
  void process();    

  class fRoot
  {
  public:
    double operator()(double x)
    {
      double delta,beta;
      delta=atan(x/g);
      beta=(x*tan((M_PI/l_box)*z))/( (M_PI/l_box)  * tan(x*z+delta));
      return (sin(x*z+delta) - pow(sin((M_PI/l_box)*z),beta));
    }
    
    double z;
    double l_box;
    double g;  
  };
  
private:
  double k2;
  fRoot fR;
  std::vector<real_t> parameters;
  
  inline real_t d0Exp(const double &x) const
  {
    return (x < parameters[3]) ? sin(parameters[0]*x + parameters[1]) : x > parameters[4] ? sin(parameters[0]*parameters[4] + parameters[1]) :  pow(sin(k2*x),parameters[2]);
  }

};

class jastrow_delta_in_trap : public jastrow<jastrow_delta_in_trap>
{
  
public:

  jastrow_delta_in_trap(real_t a_) : a(a_){};
  jastrow_delta_in_trap(const json_t & j);  
  
  inline double d0(const double & x) const {return log(x + a);}
  inline double d1(const double & x) const {return 1./(x+a);}
  inline double d2(const double & x) const {return -(1./pow( x+a,2) );}
  
  static std::string name() {return "delta_in_trap";}
private:
  real_t a;
};

class jastrow_delta_bound_state_phonons : public jastrow<jastrow_delta_bound_state_phonons>
{
public:
  jastrow_delta_bound_state_phonons(const json_t & j);

  void process(real_t);
  
  inline real_t d0(const real_t x) const {return (x<xI) ? -k*x : logA + log(pow(sin(k2*x),beta));}
  inline real_t d1(const real_t x) const {return (x<xI) ? -k : beta*k2/tan(k2*x);}
  inline real_t d2(const real_t x) const {return (x<xI) ? 0 : -k2*beta*k2/pow(sin(k2*x),2);}
  
  static std::string name() {return "delta_bound_state_phonons";}
private:
  real_t lBox;
  real_t xm;
  real_t xI;
  real_t A;
  real_t k;
  real_t k2;
  real_t beta;
  real_t g;
  real_t logA;
};

