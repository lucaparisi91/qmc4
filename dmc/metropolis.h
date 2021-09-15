#ifndef METROPOLIS_H
#define METROPOLIS_H


#include <random>
#include "traits.h"

class metropolis
{
public:
  metropolis();
  bool acceptLog(real_t ratioLog,randomGenerator_t & randGenerator);

  real_t getAcceptanceRatio() const;

  void clear();
  
  void accumulateMPI(int root);

  void recordAcceptance(bool accept ) ;
  
private:
  std::uniform_real_distribution<real_t> uniformDis;
  size_t n;
  size_t nAccepted;
};


#endif
