#include "moves/move.h"
#include <random>

class walker;

class gaussianMover : public mover
{
public:
	using walker_t=walker;

	gaussianMover(real_t sigma_);
	virtual void move(walker_t & current_walker
	,const walker_t & old_walker
	, rand_t & randGen) override; // move particle configurations
	bool isReversible() {return true;}
	real_t virtual transitionProbabilityRatio
	(const walker_t & w_new,const walker_t & w_old) {return 1.;}; /* ratio of transition probability and inverse transition probability : G(x,x')/G(x',x)
priv*/

private:
	real_t sigma;
	std::normal_distribution<real_t> distribution;
};