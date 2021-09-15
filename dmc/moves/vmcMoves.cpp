#include "vmcMoves.h"
#include "walkers.h"
#include <iostream>
#include "tools.h"

gaussianMover::gaussianMover(real_t sigma_ ) 
: sigma(sigma_),distribution(0.,1.)
{
  
}
void gaussianMover::move( 
	gaussianMover::walker_t & new_walker,
	const gaussianMover::walker_t & old_walker,
	gaussianMover::rand_t & randGenerator

	)
{
	
	auto & current_states = new_walker.getStates();
	const auto & old_states = old_walker.getStates();
	
 	for (int j=0;j<current_states.size();j++)
	{
		auto & current_state = current_states[j];
		const auto & old_state = old_states[j];

		int N=getN(current_state);
		constexpr int D= getDimensions();
		for(int i=0;i<N;i++)
			for (int d=0;d<D;d++)
		{
			current_state(i,d)=old_state(i,d) + distribution(randGenerator)*sqrt(sigma); 
		};

	}
}
