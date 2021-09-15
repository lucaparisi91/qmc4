#ifndef BRANCHING_H
#define BRANCHING_H
#include "traits.h"
class dmcWalker;
template <class T> class walkerContainer;


class branchingControl 
{
public:
  using walker_t = dmcWalker;
  using walkerContainer_t = walkerContainer<dmcWalker>;

  
  branchingControl(real_t timeStep_,real_t meanWalkers_,real_t deltaWalkers_);

  void branch(walkerContainer_t & newWalkers,const walkerContainer_t & oldWalkers,randomGenerator_t & rand);
  int nDescendants(walker_t & new_walker, const walker_t & old_walker,randomGenerator_t & rand);

  void setEnergyShift(const walkerContainer_t & walkers);
  
private:
  size_t meanWalkers;
  size_t deltaWalkers;
  size_t maxWalkers;
  real_t energyShift; // E_T , can be offset to costrain the number of walkers
  real_t timeStep;
  std::uniform_real_distribution<real_t> uniformDis;
  std::vector<int> _nDescendants;
};

#endif
