#ifndef STRUCTURE_FACTOR_H
#define STRUCTURE_FACTOR_H


#include "observables.h"


#if DIMENSIONS == 1

class structureFactor : public realVectorObservable
{
public:
  structureFactor(const json_t & j) ;
  structureFactor(int setA,int setB,const vec_t & qs,int channel=1);
  virtual const vec_t &  x() override {return _qs;}
  
  virtual std::vector<int> sets() const override {return {_setA,_setB}; };
  
  static std::string name() {return "structureFactor";}
  
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc);

  
private:
  vec_t _qs;
  int _setA;
  int _setB;
  int _channel;
};

#endif


#endif
