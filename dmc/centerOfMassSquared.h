#ifndef CENTER_OF_MASS_SQUARED_H
#define CENTER_OF_MASS_SQUARED_H
#include "observables.h"

class centerOfMassSquared : public realScalarObservable
{
public:
  centerOfMassSquared(int setA_) : setA(setA_) {};
  centerOfMassSquared(const json_t & j);
  
  static std::string name() {return "centerOfMassSquared";}
  
  virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;
private:
  difference_t diff;
  distance_t dis;
  int setA;
};

#endif
