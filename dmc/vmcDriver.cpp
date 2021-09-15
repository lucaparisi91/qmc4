#include "vmcDriver.h"
#include "walkers.h"
#include "moves/vmcMoves.h"
#include "estimators.h"
#include "timer.h"
#include "initializer.h"
#include "wavefunction/productWavefunction.h"
#include "ptools.h"
#include "tools.h"
#include <sys/stat.h>

vmcDriver::vmcDriver(vmcDriver::wavefunction_t * wave_,real_t sigma_) :
driver::driver(wave_),
 metropolisObj(),
 vmcMove(new gaussianMover(sigma_))
 {
 	
 }

void vmcDriver::run(states_t & states,size_t nBlocks)
{
	wavefunction_t & wave= getWavefunction();
	
	bool initStatus=initializer::initialize(current_walker,states,getWavefunction());
	if ( !initStatus )
	  {
	    throw invalidInput("Invalid initial conditions.");
	  }
	initializer::initialize(old_walker,states,wave);
	initializer::initialize(tmp_walker,states,wave);

	initializer::addDistances(current_walker,getEstimators());
	initializer::addDistances(old_walker,getEstimators());
	
	driver::run(nBlocks);
	
};

bool update(walker & w,productWavefunction & psi)
	{
	  w.getTableDistances().update(w.getStates());
	  w.getTableSlaters().update(w.getStates());
	  
	  if (psi.satisfyConstraints(w) )
	    {
	      w.getLogWave()=psi(w);
	      return true;
	    }
	  else
	    {
	      return false;
	    }
	};


void vmcDriver::step()
{
	auto & wave=getWavefunction();

	std::swap(old_walker,current_walker);

	vmcMove->move(current_walker,old_walker,getRandomGenerator());
	bool accept=true;
	
	// update distances and evaluates the wavefunction
	bool isValidWalker = update(current_walker,wave);
	// accept or reject the walker

	if (isValidWalker)
	  {
	    
	    accept = metropolisObj.acceptLog(2* (current_walker.getLogWave() - old_walker.getLogWave()),getRandomGenerator());
	  }
	else
	  {
	    accept=false;
	    metropolisObj.recordAcceptance(false);
	    
	  }
	if (!accept)
	{
		current_walker=old_walker;
	}
}

void vmcDriver::accumulate()
{
	auto & wave=getWavefunction();

	for( auto & est : getEstimators())
	{
		est->accumulate(current_walker,wave);
	}
}

void vmcDriver::out()
{
  auto & ests = getEstimators();
  ests.accumulateMPI(0.);
  metropolisObj.accumulateMPI(0);

  
  if (pTools::isMaster() )
    {

      std::cout << ansiColor("green") << "Block: "<< ansiColor("default")<<getCurrentBlock()<<std::endl;
      std::cout << "Acc. Ratio: " << metropolisObj.getAcceptanceRatio() << std::endl;
      
	// erase
	
	auto & energyEst = getEstimators()[0];
	std::cout << ansiColor("cyan") << "Energy: " << ansiColor("default");
	energyEst->write(std::cout);
	std::cout << std::endl;
	
	ests.dump();

    }
	ests.clear();
	metropolisObj.clear();
	confSave.dump(current_walker,getCurrentBlock());
}

void configurationsSaver::dump(const walker & w, int i)
{
  json_t j=::toJson(w.getStates() );
  
  jSonConfigurations.push_front(j);
  
  while ( jSonConfigurations.size() > nConfigurationsToSave() )
    {
      jSonConfigurations.pop_back();
    }
  
  std::ofstream f;
  int pId =pTools::rank();

  struct stat st = {0};

  if (stat(baseDir.c_str(), &st) == -1) {
    mkdir(baseDir.c_str(), 0700);
  }
  
  f.open(baseDir + "/walkers-Rank" + std::to_string(pId) + ".json");
  json_t jConfs;
  jConfs["configurations"]=json_t();
  
  for (const auto & jConf : jSonConfigurations )
    {
      jConfs["configurations"].push_back(jConf);
    }
  f << jConfs;
  f.close();
}

configurationsSaver::configurationsSaver(int nConfSave_) : baseDir("configurations"),_nConfigurationsToSave(nConfSave_)
{
  
}
