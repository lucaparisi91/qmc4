#include "branching.h"
#include "walkers.h"
#include <algorithm>
#include <iostream>

branchingControl::branchingControl(real_t timeStep_,real_t meanWalkers_,real_t deltaWalkers_) : timeStep(timeStep_),meanWalkers(meanWalkers_),deltaWalkers(deltaWalkers_),uniformDis(0.,1.),energyShift(0.),maxWalkers( (meanWalkers + deltaWalkers)*4)
{
  
}

int branchingControl::nDescendants(branchingControl::walker_t & new_walker, const branchingControl::walker_t & old_walker,randomGenerator_t & rand)
{
  return int(exp(-timeStep* ( 0.5*(new_walker.getEnergy() + old_walker.getEnergy() ) - energyShift ) ) + uniformDis(rand) );
}

void branchingControl::branch(branchingControl::walkerContainer_t & newWalkers,const branchingControl::walkerContainer_t & oldWalkers,randomGenerator_t & rand)
{
  // compute the number of descendants for each walkers
  _nDescendants.resize(newWalkers.size());

  assert(newWalkers.size() > 0);
  int newPop=0;
  
  for (int i=0;i<oldWalkers.size();i++)
    {
      _nDescendants[i]=nDescendants(newWalkers[i],oldWalkers[i],rand);
      newPop+=_nDescendants[i];
    }
  if (newPop==0)
    {
      _nDescendants[0]=1; // always keep at leat a walker
    }
  
  size_t j1=0;
  while  (j1<newWalkers.size()  )
    {
      int n=_nDescendants[j1];
      
      if (n==0) // deleter walker : swap with last walker, do not advance the index
	{
	  std::swap(newWalkers[j1],*(*(newWalkers.end() - 1 )));
	  _nDescendants[j1]=*(_nDescendants.end() - 1);
	  newWalkers.resize(newWalkers.size()-1);
	  _nDescendants.resize(_nDescendants.size()-1);
	  
	}
      else if (n>=1) // advance index and push back new copies at the end
	{

	  for(int k=0;k<n-1;k++)
	    {

	      newWalkers.push_back(newWalkers[j1]);
	      _nDescendants.push_back(1);
	    }
	  
	  _nDescendants[j1]=1;
	  j1+=1;
	}

    }

}

void branchingControl::setEnergyShift(const branchingControl::walkerContainer_t & walkers)
{
  // set shift to average energy
  energyShift=0;
  
  for (const auto & w : walkers)
    {
      energyShift+=w->getEnergy();
    }
  
  energyShift=pTools::sumAll(energyShift);
  auto population=pTools::sumAll((int)walkers.size());
  
  energyShift/=population;
  
  
  // restrict to window

  if (population > (meanWalkers + deltaWalkers))
    {
      energyShift +=1./timeStep	* log (meanWalkers*1./(meanWalkers + deltaWalkers))
	;
	
    }
  else if (population < ( meanWalkers - deltaWalkers) )
  {
    energyShift +=1./timeStep
      * log (meanWalkers*1./(meanWalkers - deltaWalkers))
	;

  }
  
}
