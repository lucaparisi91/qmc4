#include "observables.h"

class energy;
class mappedOptimizationParameter;


class optimizationObservables : public realVectorObservable
{
public:
  
  optimizationObservables(int nP,std::vector<std::vector<mappedOptimizationParameter> > parameters_ , energy & energyOb_);
  
  
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc);
  
  static std::string name() {return "optimizationObservables";}

  
private:

  
  std::vector<std::vector<mappedOptimizationParameter> > parameters;
  std::vector<real_t> parameterGradient;
  energy * energyOb;
  
};
