#ifndef MOVER_H
#define MOVER_H

#include "traits.h"
#include "qmcExceptions.h"

class walker;
class dmcWalker;
class productWavefunction;

class mover
{
public:
	/*
	Abstract class which defines move operator and log T(x',x)/T(x,x') 
	the ratio between an accepted and reject move
	*/
	mover(){}
	using rand_t = randomGenerator_t;

	virtual void move(walker & new_w,const walker & old_w,rand_t & rand_) 
	{throw missingImplementation("Move on class walker");}

	virtual real_t transitionProbabilityRatio(const walker & new_w, const walker & old_w )
	{throw missingImplementation("Transition probability on class walker");}

	virtual void move(dmcWalker & new_w,const dmcWalker & old_w,rand_t & rand_) 
	{throw missingImplementation("Move on class walker");}

	virtual real_t transitionProbabilityRatio(const dmcWalker & new_w, const dmcWalker & old_w )
	{throw missingImplementation("Transition probability on class dmcWwalker");}


private:
	
};

#endif