#include "moves/move.h"
#include <random>
#include "vmcMoves.h"
class dmcWalker;

class driftDiffusionFirstOrder : public mover
{
public:
	using walker_t = dmcWalker;
	driftDiffusionFirstOrder(real_t timeStep);

	virtual void move( 
		walker_t & new_walker,
		const walker_t & old_walker,
		rand_t & randGenerator
		);

	virtual real_t transitionProbabilityRatio(const dmcWalker & new_w, const dmcWalker & old_w );

private:
	real_t _timeStep;
	gaussianMover diffusionMove;
};