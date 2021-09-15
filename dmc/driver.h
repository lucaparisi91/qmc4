#ifndef DRIVER_H
#define DRIVER_H

#include "metropolis.h"
#include "traits.h"
#include "walkers.h"
#include "estimatorCollection.h"
#include "moves/move.h"
#include <deque>

class productWavefunction;

class estimatorBase;
class walker;

class driver
{
public:
    using estimators_t = estimatorCollection;
	using wavefunction_t = productWavefunction;
	driver(wavefunction_t * wave_);
  
	auto & getRandomGenerator() {return _rand;}

	size_t  getStepsPerBlock() const  {return _stepsPerBlock;}
        auto &  getStepsPerBlock()   {return _stepsPerBlock;}
  
        
	size_t getCurrentBlock() const {return iBlock;}

	size_t getCurrentSubStep() const  {return iSubStep;}

  size_t getCorrelationSteps() const  {return correlationSteps;}

  size_t & getCorrelationSteps()   {return correlationSteps;}
  
  
  auto & getEstimators() {return _estimators;}

  auto getCurrentStep() const  {return getCurrentBlock() * getStepsPerBlock() + getCurrentSubStep();}
  
  
	virtual void run(size_t nBlocks);
	virtual void step()=0;
	virtual void accumulate()=0;
  virtual void out()=0;
  
  virtual void isend() {};
 
  wavefunction_t & getWavefunction () {return (*_wave);}
  const wavefunction_t & getWavefunction() const {return (*_wave);}
  
private:
  wavefunction_t * _wave;
  randomGenerator_t _rand;
  size_t _stepsPerBlock;
  size_t iBlock;
  size_t iSubStep;
  estimators_t _estimators;
  size_t correlationSteps;
};






#endif
