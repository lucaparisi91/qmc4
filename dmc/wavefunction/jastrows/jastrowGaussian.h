
#include "jastrow.h"
#include "parameters.h"


class gaussianJastrow : public jastrow<gaussianJastrow>
{
	// J(x) = -alpha* x^2 
public:
  using jastrow<gaussianJastrow>::addGradientParameter;
  
  int nParameters() {return 1;}; // number of variational parameters supported

  
  gaussianJastrow(real_t alpha_) : alpha(alpha_){};
  gaussianJastrow(const json_t  & j)
  {
    alpha=j["alpha"].get<real_t>();
  }
  
  static std::string name() {return "gaussian"; }
  
  real_t d0(real_t x) const {return -alpha*x*x;}
  real_t d1(real_t x) const {return -2.*alpha*x;}
  real_t d2(real_t x) const {return -2*alpha;}
  
  void addGradientParameter(real_t x, const optimizationParameter & param, iterator_t begin, iterator_t end )
  {
    if ( param.name()=="alpha")
      {
	(*begin)-=x*x;
      }
    
  }
	
private:
  real_t alpha;
};