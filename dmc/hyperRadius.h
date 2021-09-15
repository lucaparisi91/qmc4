#ifndef HYPER_RADIUS_H
#define HYPER_RADIUS_H

#include "observables.h"



class trimerhyperRadius : public realScalarObservable
{
public:

  trimerhyperRadius(int setA_,int setB_,int setC_);
  trimerhyperRadius(const json_t & j);
  
  static std::string name() {return "trimerhyperRadius";}
  
  
  virtual std::vector<int> sets() const override {return {setA,setB,setC} ;};
  
  virtual real_t operator()(walker_t & w,wavefunction_t & psi) override;

private:
  int setA,setB,setC;
};

#endif
