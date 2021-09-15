#include "ptools.h"
#include "gtest/gtest.h"
#include <random>
#include <ctime>
#include "walkers.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include "traits.h"
#include "tools.h"
#include "wavefunction/jastrows/jastrow.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "potential.h"
#include "energy.h"
#include "initializer.h"
#include "centerOfMassSquared.h"
#include "correlationEstimator.h"
#include "wavefunction/jastrows/jastrowGaussian.h"


TEST(pToolsTest,alias_load_balancing_algorithm_test)
{
  int Np=10;
  int Nw=94;
  int nTrials=10;
  
  int k=Nw/Np;
  
  std::vector<int> permutations;
  std::vector<int> populations(Np,k);
  
  for(int i=0;i<Nw%Np;i++)
    {
      populations[i]+=1;
    }
  
    
  std::vector<int> sources;
  std::vector<int> amounts;


  std::ranlux24  randGen(time(NULL));
  std::uniform_int_distribution<int> walkerSelector(0,Np-1);
  std::uniform_int_distribution<int> amountChange(0,int(0.5*k));
  
  
  for (int i=0;i<nTrials;i++)
    {
      int n=amountChange(randGen);
      populations[walkerSelector(randGen) ]+=n;
      do
	{
	  int i=walkerSelector(randGen);
	}
      while (populations[i]<n);
      populations[ i ]-=n;
    }
  
  std::vector<int> tmpPopulations(populations);
  std::vector<std::vector<int> > destinations;
  
  pTools::determineLoadBalanceComunicationsAliasMethod(tmpPopulations,permutations,sources,destinations,amounts);

  std::vector<int> newPopulations(populations);
  
  for ( int i =0 ; i < populations.size();i++)
    {
      // increases the current population
      newPopulations[i]+= amounts[i];
      newPopulations[  sources[i]   ]-=amounts[i];
    }
  
  auto size2 = std::accumulate(newPopulations.begin(),newPopulations.end(),0);
  
  ASSERT_EQ( Nw , size2);
  for(int i=0;i<populations.size();i++)
    {
      ASSERT_NEAR (Nw/Np, newPopulations[i] ,1);
    }
  
}

TEST(pToolsTest,walkerSend)
{
  
  if ( pTools::nProcesses() >=2 )
{
  
      std::vector<int> Ns{33};
  
  real_t lBox=10000.;
  int seed=100;

  
  geometryPBC geo(lBox,lBox,lBox);
  state_t state1(Ns[0],getDimensions() );
  state_t state2(Ns[0],getDimensions());


  srand((unsigned int) seed);

  
  state1.setRandom();
  state2.setRandom();
  
  states_t states1{state1};
  states_t states2{state2};
  
  real_t alpha=1.;
  
  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction psi;
  psi.add(&wave);
  
  harmonicPotential v(geo,1.,0);
  sumPotentials pot({&v});

  energy eO(&pot);

  dmcWalker w1;
  dmcWalker w2;
  
  initializer::initialize(w1,states1,psi,eO);
  initializer::initialize(w2,states2,psi,eO);
  
  if (pTools::rank() == 1)
    {
      pTools::partialSend(w1,0,0);
    }
  else if (pTools::rank()==0)
    {
      pTools::partialRecv(&w2,1,0);
      
      auto & statew1 = w1.getStates()[0];
      auto & statew2= w2.getStates()[0];

      auto & gradw1 = w1.getGradients()[0];
      auto & gradw2 = w2.getGradients()[0];
      
      ASSERT_EQ(w1.getEnergy(),w2.getEnergy());
      ASSERT_EQ(w1.getLogWave(),w2.getLogWave());
      ASSERT_EQ(w1.getLaplacianLog(),w2.getLaplacianLog());

      for (int i=0;i<Ns[0];i++ )
	{
	  ASSERT_NEAR(statew1(i,0),statew2(i,0) , 1e-5);
	  ASSERT_NEAR(gradw1(i,0),gradw2(i,0) , 1e-5);
	}
      
      
    }
  
}
  
}

 auto getEnergies(const   pTools::walkerDistribution::walkers_t & walkers)
  {
  /* Gather all energies on process 0 for this walker distribution */
    std::vector<real_t> energies;
    for (int i=0;i<walkers.size();i++)
      {
	energies.push_back(walkers[i].getEnergy() );
    }

  std::vector<int> populations;
  populations.resize(pTools::nProcesses() , 0 );
  int n = energies.size();
  MPI_Gather( & n,1 , MPI_INT,
	      populations.data(), 1, MPI_INT,
	      0, MPI_COMM_WORLD);

  
  std::vector<real_t> totalEnergies;
  int totN = 0;
  for (int i=0;i<populations.size();i++)
    {
      totN+=populations[i];
    }
  std::vector<int> offsets;
  offsets.resize(populations.size(),0);
  
  std::partial_sum(populations.begin(),populations.end() - 1,offsets.begin() + 1 );
  totalEnergies.resize(totN);
  
  MPI_Gatherv(
	      energies.data(),
	      energies.size(),
	      MPI_DOUBLE,
	      totalEnergies.data(),
	      populations.data(),
	      offsets.data(),
	      MPI_DOUBLE,
	      0,
	      MPI_COMM_WORLD
	      );
  return totalEnergies;
  
}

TEST(pTools,walkerSend)
{
  int base_seed = 34;
  for (  int idx = 0; idx < 1000  ; idx++    )
    {
      int seed = base_seed + idx;
      
      std::vector<int> Ns{33};
  
      real_t lBox=10000.;

  
  geometryPBC geo(lBox,lBox,lBox);
  
  real_t alpha=1.;
  
  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction psi;
  
  psi.add(&wave);
  
  harmonicPotential v(geo,1.,0);
  sumPotentials pot({&v});

  energy eO(&pot);
  
  pTools::walkerDistribution wd;

  pTools::walkerDistribution::walkers_t walkers;

  int nP = pTools::nProcesses();
  
  std::ranlux24 randGen(seed + pTools::rank() );
  srand(seed + pTools::rank());
  std::uniform_int_distribution<int> disWalker(80,110);
  
  walkers.resize( disWalker(randGen) );

  
  
  
  
  for(int i=0;i<walkers.size();i++)
    {
      state_t state(Ns[0],getDimensions() );
      state.setRandom();
      
      initializer::initialize(walkers[i],{state},psi,eO);
    }

  auto oldPopulation = wd.gatherPopulations(walkers.size());
  
  auto oldEnergies = getEnergies(walkers);
  
  wd.isendReceive(walkers);
  
  wd.wait(walkers);
  
  auto newPopulation = wd.gatherPopulations(walkers.size() );
  
  auto newEnergies = getEnergies(walkers);
  
  
   if (pTools::rank()==0)
     {
       
       int oldPopulationSize=0;
       for (auto & pop : oldPopulation)
   	 {
   	   oldPopulationSize+=pop;
   	 }
       int newPopulationSize=0;
       for (auto & pop : newPopulation)
   	 {
   	   newPopulationSize+=pop;
   	 }
       
       int k = oldPopulationSize/pTools::nProcesses();
       
       for (int i=0;i<newPopulation.size();i++)
   	 {
   	   ASSERT_NEAR(newPopulation[i] ,k,1);
   	 }

       ASSERT_EQ( oldPopulationSize ,newPopulationSize);
       
       
       
       for (int i=0;i<oldEnergies.size();i++)
   	 {
   	   bool found = std::find(newEnergies.begin(),oldEnergies.end(), oldEnergies[i]) != newEnergies.end() ;
   	   ASSERT_EQ(found,1);
   	 }
       
       
     }
    }
}






TEST(pTools,iWalkerSendReceive)
{
  
   if ( pTools::nProcesses() >=2 )
     {
  
      std::vector<int> Ns{33};
  
      real_t lBox=10000.;
      int seed=100;

  
      geometryPBC geo(lBox,lBox,lBox);
      state_t state1(Ns[0],getDimensions() );
      state_t state2(Ns[0],getDimensions());


      srand((unsigned int) seed);

  
      state1.setRandom();
      state2.setRandom();
  
      states_t states1{state1};
      states_t states2{state2};
  
      real_t alpha=1.;
  
      auto J=gaussianJastrow(alpha);
      jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
      productWavefunction psi;
      psi.add(&wave);
      
      harmonicPotential v(geo,1.,0);
      sumPotentials pot({&v});

      energy eO(&pot);

      dmcWalker w1;
      dmcWalker w2;
      
      initializer::initialize(w1,states1,psi,eO);
      initializer::initialize(w2,states2,psi,eO);
      centerOfMassSquared ob(0);
      realScalarStorer storer( "cm2",new centerOfMassSquared(0) , 100 );
      centerOfMassStorer cmStorer("cm",100,{0,0});

      storer.reserve(w1);
      storer.reserve(w2);
      cmStorer.reserve(w1);
      cmStorer.reserve(w2);

      
      for(int i=0;i<355;i++)
	{
	  w1.getStates()[0].setRandom();
	  storer.store(w1,psi);
    cmStorer.store(w1,psi);
	}
      MPI_Request sendReq;
      MPI_Request recvReq;

      if ( pTools::rank() == 1 )
	{
	  pTools::ipartialSend(w1,0,0,&sendReq);
	  MPI_Status status;

	  MPI_Wait(&sendReq,&status);
	}
      else if (pTools::rank() == 0)
	{
	  pTools::ipartialRecv(&w2,1,0,&recvReq);
	  MPI_Status status;
	  MPI_Wait(&recvReq,&status);

	  auto & data2 = w2.getStorageScalarCorrelators().at("cm2");
	  auto & data1 = w1.getStorageScalarCorrelators().at("cm2");

    auto & dataCM1 = w1.getStorageScalarCorrelators().at("cm");
    auto & dataCM2 = w2.getStorageScalarCorrelators().at("cm");
    

	  for (int i=0;i<data1.size();i++)
	    {
	      EXPECT_NEAR( data1(i) ,data2(i) , 1e-5 );
	    }
	  
    for (int i=0;i<dataCM1.size();i++)
	    {
	      EXPECT_NEAR( dataCM1(i) ,dataCM1(i) , 1e-5 );
	    }
  
    ASSERT_EQ(w1.getTimeIndex(),w2.getTimeIndex());
    ASSERT_EQ(w1.getFillingStatus(),w2.getFillingStatus());
	  
	  
	}
      
     }
      
    
      
 }
