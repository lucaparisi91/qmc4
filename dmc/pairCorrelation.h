#include "observables.h"

class pairCorrelation : public realHistogramObservable
{
public:
  pairCorrelation(int setA,int setB);
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc);
  
  pairCorrelation(const json_t & j);
  
  std::vector<int> sets() {return {setA,setB};}


  static std::string name() {return "pairCorrelation";}
  
  virtual std::vector<int> sets() const override {return {setA,setB}; }
  
private:
  void setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const accumulator_t & acc);

  int setA;
  int setB;
  
    
  real_t _normalizationFactor;
  real_t max;
  real_t deltax;
};
