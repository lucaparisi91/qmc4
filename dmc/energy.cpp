#include "energy.h"
#include "wavefunction/productWavefunction.h"
#include "walkers.h"
#include <iostream>

real_t kineticEnergy::operator()(walker_t & w,wavefunction_t & psi)
	{
		real_t e=0;
		real_t ef=0;
		real_t dummy=0;
		
		psi.evaluateDerivatives(w);
		e=w.getLaplacianLog();
		
		for (const auto & grad : w.getGradients())
		{
		  real_t tmp = (grad.array() * grad.array() ).sum();
			ef+=tmp;	
		}
		
		for (const auto & grad : w.getPhaseGradients())
		{
		  real_t tmp = (grad.array() * grad.array() ).sum();
		  ef-=tmp;	
		}
		
		return -0.5*(ef + e);
	};

real_t energy::operator()(walker_t & w,wavefunction_t & psi)
	{
	  
	  auto v=(*_pot)(w);
	  return  v + kinE(w,psi);
	}; 

real_t forceEnergy::operator()(walker_t & w,wavefunction_t & psi)
	{
		real_t e=0;
		real_t ef=0;
		real_t dummy=0;
		psi.evaluateDerivatives(w);
		e+=w.getLaplacianLog();
		
		for (const auto & grad : w.getGradients())
		{
		  real_t tmp = (grad.array() * grad.array() ).sum();
			ef+=tmp;	
		}
		
		for (const auto & grad : w.getPhaseGradients())
		{
		  real_t tmp = (grad.array() * grad.array() ).sum();
		  ef+=tmp;	
		}
		
		auto v=(*_pot)(w);
		return 0.5*(ef ) + v;

		
	};

real_t energyFromWalker::operator()(energyFromWalker::walker_t & w,energyFromWalker::wavefunction_t & psi) 
{
  return w.getEnergy();
}
