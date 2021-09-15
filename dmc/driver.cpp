#include "driver.h"
#include "wavefunction/productWavefunction.h"
#include "estimators.h"
#include "walkers.h"
#include "initializer.h"
#include <iostream>
#include "tools.h"
#include "moves/vmcMoves.h"
#include "ptools.h"
#include  <filesystem>

#include "timer.h"
#include <sys/stat.h>


driver::driver(driver::wavefunction_t * wave) : _wave(wave), 
						iBlock(0),iSubStep(0),_stepsPerBlock(0),correlationSteps(0)
{
  getTimers().add("totalRunningTime");
  getTimers().add("outputTime");
}

void driver::run(size_t nBlocks)
{
  getTimers().get("totalRunningTime").start();
	for (iBlock=0;iBlock<nBlocks;iBlock++)
	{
		for (iSubStep=0;iSubStep<_stepsPerBlock;iSubStep++)
		{
		  for (int iCorrelationStep=0;iCorrelationStep<correlationSteps;iCorrelationStep++)
		    {
		      
		      step();
		      isend();
		    }
		  step();
		  accumulate();
		  isend();
		};
		getTimers().get("outputTime").start();
		out();
		getTimers().get("outputTime").stop();
				
	}
  getTimers().get("totalRunningTime").stop();

  if (pTools::rank() == 0)
    {
      auto totalRunningTime = getTimers().get("totalRunningTime").milliseconds();
      auto it = getTimers().begin();
      while( it != getTimers().end() )
	{
	  auto elapsed = (it->second).milliseconds();
	  std::cout<<it->first << " :  "<<  elapsed <<" ms ( "<< elapsed*1./totalRunningTime * 100 <<" %  )"<< std::endl;
	  it++;
	}
    }
  
};

