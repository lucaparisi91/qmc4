#include "parameters.h"
#include "wavefunction/jastrows/jastrow.h"
#include <iostream>
#include <unsupported/Eigen/CXX11/Tensor>
#include "geometry.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "potential.h"
#include "initializer.h"
#include "walkers.h"
#include "energy.h"
#include "estimators.h"
#include "dmcDriver.h"
#include "vmcDriver.h"
#include "branching.h"
#include <string>
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "factory.h"
#include "wavefunction/jastrows/jastrowSquareWell.h"
#include "wavefunction/jastrows/jastrow_delta.h"
#include "wavefunction/jastrows/jastrowPoschTeller.h"
#include "wavefunction/jastrows/jastrowSpline.h"
#include "wavefunction/jastrows/jastrowGaussian.h"


#include "ptools.h"
#include "pairCorrelation.h"
#include "centerOfMassSquared.h"
#include "hyperRadius.h"
#include "wavefunction/jastrows/jastrowDipolar.h"
#include "structureFactor.h"
#include "trimerCorrelation.h"
#include "wavefunction/jastrows/jastrowGaussGauss.h"
#include "potentialFunctors.h"

#include "wavefunction/jastrows/jastrowHardSphereGauss.h"


bool check_n_particles(const states_t & states,const std::vector<int> & Ns)
{
  bool pass=true;
  for (int i=0;i<Ns.size();i++)
    {
      pass=pass and ( getN(states[i]) == Ns[i] );
    }
  return pass;
}

bool containSuperfluidFractionMeasurements( const json_t & j)
{
  if (
      j["measurements"].find_if(
				[](const json_t & jm) { return jm["kind"]=="superfluidFraction" ;}
				)
      != j["measurements"].end()
      )
	   {
	     return true;
	   }
  else
    {
      return false;
    }
}




auto generateRandomInitialConditions(const std::vector<int> Ns, const std::vector<real_t> & lBoxInitialCondition)
{
  states_t states;
  const int D = lBoxInitialCondition.size();
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , D);
      particleData.setRandom();

      for (int i=0;i<getN(particleData);i++)
	{
	  for (int d=0;d<getDimensions();d++)
	    {
	      particleData(i,d)*=lBoxInitialCondition[d]*0.5;
	    }
	}
      //particleData*=lBox[0]/2.;
      states.push_back(particleData);
    }
  return states;
}

std::string readToString(istream & inStream)
{
  std::string jSonString = "";
  std::string line;

  while (std::getline(inStream, line))
    {
      jSonString+=line;
    }
  return jSonString;

}


int main(int argc, char** argv)
{

  /*
    Reads input from standrard input and broadcast to all other mpi processes
*/
  pTools::init(argc,argv);
  std::string jSonString="";  

    if ( pTools::isMaster())
      {

        if (argc==1)
        {
          jSonString=readToString(std::cin); 
        }
        else
        {
          ifstream fd(argv[1]);
          jSonString = readToString(fd);
        }
        
    
      }
    
  

  pTools::broadcast(&jSonString,0);
  
  std::istringstream ss(jSonString);
  
  json_t j;
  
  ss >> j;

  
  

  if (pTools::isMaster())
     {
       std::cout << "MPI processes: " <<  pTools::nProcesses() << std::endl;
     }
  
  
  std::vector<real_t> lBox;
  std::vector<real_t> lBoxInitialCondition;
  lBox=j["lBox"].get<decltype(lBox)>();
  std::vector<int> Ns;
  Ns=j["N"].get<decltype(Ns)>();
  int D=lBox.size();
  
  if (j.find("initialConditionGenerator") != j.end() )
    {
      const auto & confGenJ = j["initialConditionGenerator"];
      if (confGenJ.find("lBox") != confGenJ.end() )
	{
	  lBoxInitialCondition=confGenJ["lBox"].get<decltype(lBoxInitialCondition)>() ;
	}
      else
	{
	  lBoxInitialCondition=lBox;
	}
	
    }
  else
    {
      lBoxInitialCondition = lBox;
    }
  if ( ( D != getDimensions() )or (lBoxInitialCondition.size() != getDimensions() ) )
    {
      throw invalidInput("Input file implies dimensionality different from " + std::to_string(getDimensions() ) );
      
    }
  
#if DIMENSIONS == 3
  geometryPBC geo( lBox[0], lBox[1], lBox[2]);
#endif
#if DIMENSIONS == 1
  geometryPBC geo( lBox[0] );
#endif

  
  
  getFactory().registerJastrow< gaussianJastrow >();
  getFactory().registerJastrow< jastrowGaussGauss >();
  getFactory().registerJastrow< jastrowSquareWell >();
  getFactory().registerJastrow< jastrowSquareWellBoundState >();
  getFactory().registerJastrow< jastrowPoschTeller >();
  getFactory().registerJastrow< jastrowHardSphereGauss >();
  getFactory().registerJastrow< jastrowHardSphere>();
  getFactory().registerJastrow< jastrowSpline>();

  
  #if DIMENSIONS == 3
  getFactory().registerOrbital<sinOrbital>();
  getFactory().registerOrbital<planeWave>();
  #endif
  getFactory().registerObservable<pairCorrelation>();
  getFactory().registerObservable<centerOfMassSquared>();
  getFactory().registerObservable<trimerhyperRadius>();
  getFactory().registerObservable<trimerCorrelation>();
  getFactory().registerObservable<trimerCorrelationUnDis>();
  getFactory().registerObservable<trimerCorrelationDis>();

#if DIMENSIONS == 1
  getFactory().registerObservable<structureFactor>();
#endif
  
  getFactory().registerJastrow<jastrow_delta_phonons>();
  getFactory().registerJastrow<jastrow_delta_in_trap>();
  getFactory().registerJastrow<jastrow_delta_bound_state_phonons>();
  getFactory().registerJastrow<jastrowDipolarRep>();
  
  auto waves = getFactory().createWavefunctions( j["wavefunctions"],geo);

  for (int i=0;i<waves.size();i++)
    {
      std::string filename = "wave" + std::to_string(i) + ".dat";
      ofstream f;
      f.open(filename);
      f << waves[i]->print() ;
      f.close();
      
    }
  
  productWavefunction psi(waves);
  
  /* Potentials
     First register implemented concrete potential classes and then create the local potential as sum of individual potentials
*/

  getFactory().registerPotential<harmonicPotential>();
  getFactory().registerPotential<squareWellPotential2b>();
  
  getFactory().registerPotential<dipolarPotential2b>();
  getFactory().registerPotential<poschTellerPotential2b>();

  getFactory().registerPotentialFunctor<squareWell>();
  getFactory().registerPotentialFunctor<potentialBarrier>();
  
  auto potentials = getFactory().createPotentials(j["potentials"],geo);

  
  sumPotentials pot(potentials);
  
  // generate a random initial conditions which satisfy the contraints imposed by the wavefunction
  
  states_t states ;

  {
    bool isSatisfied=false;
    while (! isSatisfied)
    {
      walker tmpWalker;
      states= generateRandomInitialConditions(Ns,lBoxInitialCondition);
      initializer::initialize(tmpWalker,states,psi);
      
      isSatisfied=psi.satisfyConstraints(tmpWalker);
    }

  
  }

  
  auto eO=new  energy(&pot);
  auto efO= new forceEnergy(&pot);
  
  auto m = new realScalarEstimator("energy",eO);
  auto m2= new realScalarEstimator("forceEnergy",efO);
  
  std::string method = j["method"].get<std::string>();
  

  if ( method == "vmc" and containSuperfluidFractionMeasurements(j) )
    {
      throw invalidInput("Cannot have vmc and superfluid fraction measurements");
    }

  real_t timeStep = j["timeStep"].get<real_t>();
  size_t stepsPerBlock = j["stepsPerBlock"].get<int>();
  size_t nBlocks = j["nBlocks"].get<int>();  
  int seed= j["seed"].get<int>();
  size_t correlationSteps=j["correlationSteps"].get<int>();

  
  std::vector<states_t> configurations;
  if (j.find("configurations") != j.end())
    {
      configurations=readStatesFromDirectory(j["configurations"].get<std::string>());
      if (
	  (configurations.size() == 0)  and
	  ( (j.find("initialConfigurations") != j.end()) )

	  )
	{
	  configurations=readStatesFromDirectory(j["initialConfigurations"].get<std::string>() );
	}
      
      if(  (pTools::rank() == 0 ) and ( configurations.size() == 0) )
	{
	  std::cout << ansiColor("yellow") << "WARNING: Failed to load configurations from directory. Falling back to randomly generated configurations." << ansiColor("default") << std::endl;
	}
      
    }

  auto ests = getFactory().createEstimators(j["measurements"]);
  auto storers = getFactory().createStorers(j["measurements"]);

  
  
  if (method == "vmc")
    {
      vmcDriver vmcO(&psi,timeStep);
      vmcO.getStepsPerBlock()=stepsPerBlock;
      vmcO.getCorrelationSteps()=correlationSteps;
      

      vmcO.getEstimators().push_back(m);

      	  if (
	      j["measurements"].find_if(
					[](const json_t & jm) { return jm["kind"]=="forceEnergy" ;}
				)
      != j["measurements"].end()
	      )      
      // if ( find(j["measurements"],"kind","forceEnergy") != j["measurements"].end() )
      	{
      	  vmcO.getEstimators().push_back(m2);
      	}
      
      for (auto & est : ests)
	{
	  vmcO.getEstimators().push_back(est); 
	}
      
      
      vmcO.getRandomGenerator().seed(seed + pTools::rank() );
      states_t * initialConfiguration = &states;
      if (configurations.size() > 0 )
	{
	  initialConfiguration = &(configurations[0]);
	}
  
      if (! check_n_particles(*initialConfiguration,Ns) )
	{
	  throw invalidInput("Initial configuration does not math the numper of particles defined in the input file");
	}

      int nConfigurationsToSaveVMC = json_t::load(j["nConfigurationsToSaveVMC"],1);





      vmcO.getConfigurationsSaver().nConfigurationsToSave() = nConfigurationsToSaveVMC;
      vmcO.run(*initialConfiguration,nBlocks); 
    }
  
  else if ( method == "dmc" or method == "svmc")
    {
      size_t nW=j["walkers"].get<int>();
      size_t delta_walkers=size_t(nW*0.1);
      
      
      
      if (j.find("delta_walkers")!= j.end() )
	{
	  delta_walkers=j["delta_walkers"].get<int>();
	}
      dmcDriver dmcO(&psi,&pot,timeStep,nW,delta_walkers);
      dmcO.getStepsPerBlock()=stepsPerBlock;
      dmcO.getCorrelationSteps()=correlationSteps;


      std::vector<states_t> dmcStates(configurations);      

      
      
      if ( dmcStates.size() == 0 )
	{
	  for(int i=0;i<nW/pTools::nProcesses();i++)
	    {
	      dmcStates.push_back(states);
	    }
      
	}

      
      dmcO.getRandomGenerator().seed(seed + pTools::rank() );
      
      if (method == "svmc")
	{
	  dmcO.disableBranching();
	  
	  if (
	      j["measurements"].find_if(
					[](const json_t & jm) { return jm["kind"]=="forceEnergy" ;}
				)
      != j["measurements"].end()
	      )

	{
	  dmcO.getEstimators().push_back(m2);
	}
		
	}

      
      for (auto & st : storers)
	{
	  dmcO.getEstimators().push_back(st);
	}
      
      
      for (auto & est : ests)
	{
	  dmcO.getEstimators().push_back(est);
	}
      

      dmcO.run(dmcStates,nBlocks);
    }

  pTools::finalize();
  
}
