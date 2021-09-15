#include "walkers.h"
#include "dmcDriver.h"
#include "wavefunction/productWavefunction.h"
#include "initializer.h"
#include "tools.h"
#include "moves/dmcMoves.h"
#include "estimators.h"
#include "branching.h"
#include "timer.h"

bool noMetropolisPolicy::accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t &generator)
{
	return true;
}

real_t noMetropolisPolicy::getAcceptanceRatio() const
{
  return 1.;
}


bool metropolisPolicy::accept(
		mover & move,
		const dmcWalker & new_walker,
		const dmcWalker & old_walker,
		wavefunction_t & wavefunction,
		randomGenerator_t &generator)
{
  real_t transitionRatio = move.transitionProbabilityRatio(new_walker,old_walker);

  return metropolisSampler.acceptLog(2*(new_walker.getLogWave() - old_walker.getLogWave() )  + transitionRatio ,generator);
  
}

real_t metropolisPolicy::getAcceptanceRatio() const
{
  return metropolisSampler.getAcceptanceRatio();
}

void metropolisPolicy::clear() 
{
  return metropolisSampler.clear();
}

void metropolisPolicy::accumulateMPI(int root) 
{
  return metropolisSampler.accumulateMPI(root);
}


void dmcDriver::update(dmcWalker & wNew, dmcWalker & wOld)
{
  dmcMover->move(wNew,wOld,getRandomGenerator());
  bool accepted;
  bool isValid=updateForceGradientEnergy(wNew, getWavefunction(),energyOb);

  
  if (isValid)
    {
      accepted=accepter->accept(*dmcMover,wNew,wOld,getWavefunction() ,getRandomGenerator());
    }
  else
    {
      accepted=false;
      accepter->recordAcceptance(false);
    }
	  
  if (!accepted)
    {
      wNew=wOld;
    }
  else
    {
      std::swap(wOld.getStorageScalarCorrelators() ,wNew.getStorageScalarCorrelators() );
      std::swap(wOld.getTimeIndex() ,wNew.getTimeIndex() );
      std::swap(wOld.getFillingStatus() ,wNew.getFillingStatus() );
    }
  
}


 void dmcDriver::step()
{
  bool allKilled=false;
  
  std::swap(current_walkers,old_walkers);

  // if ( old_walkers.size() == 0 )
  //   {
  //     allKilled=true;
  //     std::cout << ansiColor("red") << "All old walkers were killed on" << pTools::rank() << " at step "<< getCurrentStep() <<ansiColor("default")<< std::endl << std::flush;
	   
  //   }

	
  current_walkers.resize(old_walkers.size(),old_walkers[std::max(old_walkers.size()-1,(size_t)0)]);

	for (int i=0;i<old_walkers.size();i++)
	{
	  update(current_walkers[i],old_walkers[i]);
	}
	

	auto oldSize=old_walkers.size();
	
	getTimers().get("waitWalkers").start();
	walkerLoadBalancer->wait(old_walkers);	
	getTimers().get("waitWalkers").stop();
	
	// if (allKilled)
	//   {
	//     std::cout << ansiColor("yellow") << "old walkers resize on" << pTools::rank()<<" to " << old_walkers.size() <<ansiColor("default")<< std::endl << std::flush;

	//   }

	
	current_walkers.resize(old_walkers.size(),old_walkers[std::max(old_walkers.size()-1,(size_t)0)]);


	// if (allKilled)
	//   {
	//     //std::cout << ansiColor("yellow") << "current walkers resize on" << pTools::rank()<<" to " << current_walkers.size() <<ansiColor("default")<< std::endl << std::flush;

	//   }

	
	
	for (int i=oldSize;i<old_walkers.size();i++)
	{
	  update(current_walkers[i],old_walkers[i]);
	}
	
	
	if (performBranching)
	  {
	    brancher->branch(current_walkers,old_walkers,getRandomGenerator());
	    getTimers().get("shiftEnergy").start();
	    brancher->setEnergyShift(current_walkers);
	    getTimers().get("shiftEnergy").stop();    
	  }
	
}


void dmcDriver::isend()
{
  
  getTimers().get("sendWalkers").start();
  walkerLoadBalancer->isendReceive(current_walkers);
  getTimers().get("sendWalkers").stop();
	
}


void dmcDriver::run( const std::vector<states_t> &states , size_t nBlocks )
{
  
  bool status=initializer::initialize(current_walkers,states,getWavefunction(),energyOb);


  if (!status)
    {
      throw invalidInput("Inalid initial configurations.");
      
    }
  
  
  initializer::initialize(old_walkers,states,getWavefunction(),energyOb);
  
  for (auto & w : current_walkers)
    {
      getEstimators().reserve(*w);
      initializer::addDistances(*w,getEstimators());

    }
  
  for (auto & w : old_walkers)
    {
      getEstimators().reserve(*w);
      initializer::addDistances(*w,getEstimators());

    }
  
  
  //initializer::initialize(walkerLoadBalancer,getWavefunction(),energyOb);
  brancher->setEnergyShift(current_walkers);
  walkerLoadBalancer->isendReceive(current_walkers);
  
  driver::run(nBlocks);
}

void dmcDriver::out()
{
  walkerLoadBalancer->wait(current_walkers);

  auto & ests = getEstimators();
  ests.accumulateMPI(0.);
  accepter->accumulateMPI(0.);
  auto walkers_size = pTools::sum(current_walkers.size() , 0) ;

  if (pTools::isMaster() )
    {
      std::cout << ansiColor("green") << "Block: "<< ansiColor("default")<<getCurrentBlock()<<std::endl;
	
      auto & energyEst  = getEstimators()[0];
      std::cout << "Acc. Ratio: " << accepter->getAcceptanceRatio() << std::endl;
      std::cout << "Step: " << (getCurrentBlock()+1)*getStepsPerBlock()  << std::endl;

      std::cout << ansiColor("cyan") << "Energy: " << ansiColor("default");
      std::cout << std::scientific;
      energyEst->write(std::cout);
      std::cout << std::endl<<std::defaultfloat;

      std::cout << ansiColor("cyan") << "Curr. Walkers: " << ansiColor("default") ;
      std::cout << walkers_size << std::endl;
      
	
      
      ests.dump();
    }
  ests.clear();
  current_walkers.dump(getCurrentBlock());
  accepter->clear();

}

dmcDriver::dmcDriver(dmcDriver::wavefunction_t * wave, potential_t * pot,real_t timeStep,size_t nWalkers,size_t delta_walkers) :
  driver::driver(wave), energyOb(pot),
  dmcMover( new driftDiffusionFirstOrder(timeStep)),
  energyEst(new realScalarEstimator("energy",new energyFromWalker) ),
  accepter(new metropolisPolicy),
  brancher(
	   new branchingControl(timeStep,nWalkers,delta_walkers)
	   ),
  walkerLoadBalancer(new pTools::walkerDistribution ),
  performBranching(true)
  
{
  
  getEstimators().push_back( energyEst );
  getTimers().add("shiftEnergy");
  getTimers().add("sendWalkers");
  getTimers().add("waitWalkers");
  
}

void dmcDriver::disableBranching()
{
  performBranching=false;
}

void dmcDriver::accumulate()
{
  auto & wave=getWavefunction();
	  for (auto & current_walker : current_walkers)
	    {
	      getEstimators().accumulate(*current_walker,wave);
	    }
}
	


