#include "dmcMoves.h"
#include "walkers.h"
#include <cmath>
#include <iostream>

driftDiffusionFirstOrder::driftDiffusionFirstOrder(real_t timeStep) : 
_timeStep(timeStep),diffusionMove(timeStep)
{
  
}

void driftDiffusionFirstOrder::move(
	driftDiffusionFirstOrder::walker_t & new_walker ,
	const driftDiffusionFirstOrder::walker_t & old_walker,
	driftDiffusionFirstOrder::rand_t & randGenerator)
{
	// drift step
	auto & new_positions = new_walker.getStates();
	const auto & old_positions = old_walker.getStates();

	const auto & old_gradients = old_walker.getGradients();
	
	for(int i=0;i<new_positions.size();i++)
	{
	  new_positions[i]= old_positions[i] + old_gradients[i]*_timeStep; // drift
	}

       diffusionMove.move(new_walker,new_walker,randGenerator); // diffusion
       
}

real_t driftDiffusionFirstOrder::transitionProbabilityRatio(const dmcWalker & new_w, const dmcWalker & old_w )
{
	const auto & new_positions = new_w.getStates();
	const auto & new_driftForces = new_w.getGradients();

	const auto & old_positions = old_w.getStates();
	const auto & old_driftForces = old_w.getGradients();
	real_t q=0;

	

	for (int i=0;i< new_positions.size();i++)
	{
	  auto a=  (
		    _timeStep*1./2 *(old_driftForces[i] - new_driftForces[i]) 
			- new_positions[i] + old_positions[i]
		    ) ;
	  auto b=( old_driftForces[i] + new_driftForces[i] ) ;

	  real_t res= (a.array() * b.array() ).sum();
	  q+=res;
	};
	
	return q;

}
