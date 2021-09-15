#include "gtest/gtest.h"
#include <vector>
#include "dmcDriver.h"
#include "initializer.h"
#include "geometry.h"
#include "wavefunction/jastrows/jastrow.h"
#include "wavefunction/productWavefunction.h"
#include "branching.h"
#include <algorithm>
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "estimators.h"
#include "wavefunction/jastrows/jastrowGaussian.h"


TEST(dmcTest,branching_lower)
{
	int N=100;
	int D=getDimensions();
	int nW=100;	
	
 	geometryPBC geo( 10., 10., 10.);

	
 	tableDistances tab(geo);
 	real_t alpha=1;
 	auto J=gaussianJastrow(alpha);
	
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
 	productWavefunction psi({&wave});
	
 	harmonicPotential v(geo,1.,0);
	sumPotentials pot({&v});

 	energy eO(&pot);

 	walkerContainer<dmcWalker> walkers;
	walkerContainer<dmcWalker> newWalkers;

	std::vector< states_t> walkerStates;
	
	for (int i=0;i<nW;i++)
	  {
	    state_t statetmp(N,D);
	    statetmp.setRandom();
	    walkerStates.push_back ( {statetmp}); 
	  }

	initializer::initialize(walkers,walkerStates,psi,eO);
	initializer::initialize(newWalkers,walkerStates,psi,eO);

	branchingControl bc(1E-3,40,10);
	bc.setEnergyShift(walkers);
	std::vector<int> nDescendants;
	
	size_t seed= 10;
	randomGenerator_t rand(seed);
	
	nDescendants.resize(walkers.size());

	for (int i=0;i<walkers.size();i++)
	    {
	      nDescendants[i]=bc.nDescendants(newWalkers[i],walkers[i],rand);
	    }
	
	int newPopulation=std::accumulate(nDescendants.begin(),nDescendants.end(), decltype(nDescendants)::value_type(0)) ;
	
	randomGenerator_t rand2(seed);
	bc.branch(newWalkers,walkers,rand2);
	
	ASSERT_EQ( newPopulation , newWalkers.size());
	
	for (int i=0;i<walkers.size();i++)
	  {
	    real_t e=walkers[i].getEnergy();
	    
	    int n=std::count_if(newWalkers.begin(),newWalkers.end(),[&](std::unique_ptr<dmcWalker> & w){return w->getEnergy() == e;});
	    
	    ASSERT_EQ( nDescendants[i] , n );
	  }
	
 }

TEST(dmcTest,branching_upper)
{
	int N=100;
	int D=getDimensions();
	int nW=100;	
	
 	geometryPBC geo( 10., 10., 10.);
	
 	tableDistances tab(geo);
 	real_t alpha=1;
 	auto J=gaussianJastrow(alpha);
	
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);

 	productWavefunction psi({&wave});
	
	
 	harmonicPotential v(geo,1.,0);
	sumPotentials pot({&v});

 	energy eO(&pot);

 	walkerContainer<dmcWalker> walkers;
	walkerContainer<dmcWalker> newWalkers;

	
	std::vector< states_t> walkerStates;

	for (int i=0;i<nW;i++)
	  {
	    state_t statetmp(N,D);
	    statetmp.setRandom();
	    walkerStates.push_back ( {statetmp}); 
	  }

	
	
	initializer::initialize(walkers,walkerStates,psi,eO);
	initializer::initialize(newWalkers,walkerStates,psi,eO);

	branchingControl bc(1E-3,800,10);
	bc.setEnergyShift(walkers);
	std::vector<int> nDescendants;
	
	size_t seed= 10;
	randomGenerator_t rand(seed);
	
	nDescendants.resize(walkers.size());

	for (int i=0;i<walkers.size();i++)
	    {
	      nDescendants[i]=bc.nDescendants(newWalkers[i],walkers[i],rand);
	    }
	
	int newPopulation=std::accumulate(nDescendants.begin(),nDescendants.end(), decltype(nDescendants)::value_type(0)) ;
	
	randomGenerator_t rand2(seed);
	bc.branch(newWalkers,walkers,rand2);
	
	ASSERT_EQ( newPopulation , newWalkers.size());
	//EXPECT_GT(newPopulation,nW);
	for (int i=0;i<walkers.size();i++)
	  {
	    real_t e=walkers[i].getEnergy();

	    int n=std::count_if(newWalkers.begin(),newWalkers.end(),[&](std::unique_ptr<dmcWalker> & w){return w->getEnergy() == e;});
	    
	    ASSERT_EQ( nDescendants[i] , n );
	  }	
 }
