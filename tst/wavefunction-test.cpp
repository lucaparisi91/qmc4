#include "wavefunction/productWavefunction.h"
#include "gtest/gtest.h"
#include "tableDistances.h"
#include "walkers.h"
#include "initializer.h"
#include "potential.h"
#include "energy.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "wavefunction/jastrowWavefunctionTwoBody.h"

#include "wavefunction/slaterDeterminant.h"
#include "orbitals.h"
#include "estimators.h"
#include "wavefunction/jastrows/jastrowSquareWell.h"

#include "centerOfMassSquared.h"
#include "correlationEstimator.h"
#include "wavefunction/jastrows/jastrowGaussian.h"


TEST(wavefunctionTest,oneBody)
{
	int N=100;
	int D= getDimensions();
	
 	state_t particleData(N , D);
 	state_t gradient(N , D);
 	real_t alpha=1.;
 	
 	particleData.setRandom();
	
 	geometryPBC geo( 10., 10., 10.);
	

 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
 	productWavefunction waveT;
 	waveT.add(&wave);
	
 	real_t e=0 , ef=0 , waveValue =0;
	
 	states_t states {particleData};
  
	walker w;
	initializer::initialize(w,states,waveT);
 	waveT.evaluateDerivatives( w);
	e=w.getLaplacianLog();
	
 	EXPECT_EQ(e,-2*alpha*N*D);

 	auto psi_value = waveT(w);
	
 	real_t sum2 = (states[0].array()*states[0].array()).sum();
 	
 	EXPECT_NEAR(psi_value,sum2*(-alpha),1e-5);

}

TEST(wavefunctionTest,tableDistances)
{
	int N=100;
	int D=getDimensions();
 	state_t particleData(N , D);
 	state_t gradient(N , D);
	

 	particleData.setRandom();

 	geometryPBC geo( 10., 10., 10.);

 	states_t states {particleData,2*particleData};
 	tableDistances tab(geo);

 	tab.add(0);

 	tab.add(0,0);
 	tab.add(0,1);
 	tab.update(states);

 	const auto & differences = tab.differences(0,0);

 	int k=0;
 	for (int i=0;i<N;i++)
 	{
 		for(int j=0;j<i;j++)
 		{
 			for(int d=0;d<D;d++)
 			{
 				auto diff1=geo.difference(particleData(i,d) - particleData(j,d),d ) ;
 				auto diff2 = differences(k,d);

 				EXPECT_NEAR(diff1, diff2, 1e-4);	

 			}
 				k++;
 		}

 	}

}

TEST(wavefunctionTest,2b)
{
#if DIMENSIONS == 3
	int N=2;
	int D=getDimensions();
 	state_t particleData(N , D);
 	state_t gradient(N , D);
	real_t lBox=100;
	
 	particleData.setRandom();
       
	
 	geometryPBC geo( lBox,lBox,lBox);
	
 	states_t states {particleData};
	real_t aInverse=0;
	real_t V0=2.4674011002723395 ;
	real_t R0=1.0;
	squareWellPotential2b v(geo, V0 , R0, 0., 0. ) ;

	potential_t pot({&v});
	real_t Rm=40.0;
	real_t alpha=0.5;
	
        jastrowSquareWell J(V0,R0,Rm,alpha,aInverse,lBox);

	jastrowTwoBodyWavefunctionUndistinguishable<jastrowSquareWell> wave(J,geo);
	productWavefunction psi({&wave});
	dmcWalker w;

	

	energy eO(&pot);
	states[0]*=Rm/(2.*sqrt(3));
	
	initializer::initialize(w,states,psi,eO);

	double sum=0;
	const auto & dis = w.getTableDistances().distances(0,0);
	for ( int i=0;i<dis.size();i++  )
	  {
	    if (dis(i) <= R0 ) sum+=-V0;
	  }

	ASSERT_EQ(dis.size(),( N*(N-1) )/2);
	auto sum2 = v(w);
	
	EXPECT_NEAR(sum,sum2,1e-5);
	EXPECT_NEAR(eO(w,psi),0,1e-7 ) ;

#endif
	
 }



TEST(fermions,slaterWavefunctionEnergy)
{
  #if DIMENSIONS == 3
  std::vector<int> Ns{33};

  real_t lBox=33.;
  
  geometryPBC geo(lBox,lBox,lBox);
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , getDimensions());
      particleData.setRandom();
      particleData=particleData*lBox - lBox/2.;
      states.push_back(particleData);
    }
  
  real_t alpha=1.;
  
  orbitalSet<sinOrbital> sineCosBasis;
  
  fillFermiSea(sineCosBasis.getOrbitals(),Ns[0],lBox); // fills a fermi see with the given orbital based
  
  slaterDeterminantWavefunction<decltype(sineCosBasis)> wave(sineCosBasis,geo,0);

  productWavefunction psi({&wave}) ;
  dmcWalker w;
  
  emptyPotential v(geo);
  sumPotentials pot({&v});
  
  energy  eO(&pot);
  forceEnergy efO(&pot);
  
  
  initializer::initialize(w,states,psi,eO);
   

  EXPECT_NEAR( w.getEnergy() , sineCosBasis.energy() , 1E-5 );

  #endif
}


TEST(fermions,slaterWavefunctionComplexEnergy)
{
  #if DIMENSIONS == 3
  std::vector<int> Ns{33};

  real_t lBox=33.;
  
  geometryPBC geo(lBox,lBox,lBox);
  states_t states;
  
  for (int i=0;i<Ns.size();i++)
    {
      state_t particleData(Ns[i] , getDimensions());
      particleData.setRandom();
      particleData=particleData*lBox - lBox/2.;
      states.push_back(particleData);
    }
  
  real_t alpha=1.;
  
  orbitalSet<planeWave> basis;
  real_t teta=3.1;
  fillFermiSea(basis.getOrbitals(),Ns[0],lBox,teta); // fills a fermi see with the given orbital based

  
  slaterDeterminantWavefunctionComplex<decltype(basis)> wave(basis,geo,0);
  
  productWavefunction psi({&wave}) ;
  dmcWalker w;
  
  emptyPotential v(geo);
  sumPotentials pot({&v});
  
  energy eO(&pot);
  forceEnergy efO(&pot);
  
  
  initializer::initialize(w,states,psi,eO);
  
  EXPECT_NEAR( w.getEnergy() , basis.energy() , 1E-5 );

  #endif
  
}


TEST(measurementsTest,forwardWalkingCm2)
{
	
  int N=100;
  int D= getDimensions();
  state_t particleData(N , D);
  state_t gradient(N , D);
  real_t alpha=1.;
  
  particleData.setRandom();

  geometryPBC geo( 10., 10., 10.);


  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction waveT;
  waveT.add(&wave);
  
  real_t e=0 , ef=0 , waveValue =0;
	
  states_t states {particleData};
	
  walker w;
  initializer::initialize(w,states,waveT);
  centerOfMassSquared ob(0);
  realScalarStorer storer( "cm2",new centerOfMassSquared(0) , 100 );
  std::vector<real_t> cm2s;
  
  storer.reserve(w);
  for(int i=0;i<200;i++)
    {
      state_t particleData(N , D);
      particleData.setRandom();
      initializer::initialize(w,{particleData},waveT);
      storer.store(w,waveT);
      
      cm2s.push_back(ob(w,waveT) );
      
    }

  auto & storedCorrelators = w.getStorageScalarCorrelators()["cm2"];

  for (int i=0; i< std::min((int)storedCorrelators.rows() ,(int) cm2s.size() ) ; i++)
    {
      auto j = wrapIndex(w.getTimeIndex().at("cm2") -1 - i,storedCorrelators.size());
      EXPECT_NEAR(storedCorrelators(j) , cm2s[cm2s.size()-i-1] ,1e-5 )  ;
      
    }  


  walker w2;

  w2=w;

  auto & storedCorrelators1 = w.getStorageScalarCorrelators()["cm2"];
  auto & storedCorrelators2 = w2.getStorageScalarCorrelators()["cm2"];
  
  for (int i=0; i< storedCorrelators.rows() ; i++)
    {
      EXPECT_NEAR(storedCorrelators1(i) , storedCorrelators2(i) ,1e-7 )  ;

    }  
  EXPECT_EQ(w.getTimeIndex() , w2.getTimeIndex() )  ;
  
}



TEST(measurementsTest, superFluidFraction)
{
  
  int N=10;
  constexpr int D= getDimensions();
  state_t particleData(N , D);
  state_t gradient(N , D);
  real_t alpha=1.;
  
  particleData.setRandom();

  geometryPBC geo( 10., 10., 10.);


  auto J=gaussianJastrow(alpha);
  jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
  productWavefunction waveT;
  waveT.add(&wave);
  
  real_t e=0 , ef=0 , waveValue =0;
	
  states_t states {particleData};
	
  walker w;
  initializer::initialize(w,states,waveT);
  std::string storerLabel="cmSuperfluidFraction";
  
  const int nStores = 234;
  const int recordSteps= 10;
  
  Eigen::Array<real_t,Eigen::Dynamic,D> positions(nStores,D);
  
  positions.setConstant(0);
  centerOfMassStorer storer(storerLabel,recordSteps,{0});
  storer.reserve(w);

  superfluidFractionEstimator sP("rho",storerLabel,recordSteps,0,0);
  
  for(int i=0;i<nStores;i++)
    {
      state_t particleData(N , D);
      particleData.setRandom();
      initializer::initialize(w,{particleData},waveT);
      
      for (int id=0;id<D;id++)
	{
	  real_t xcmd=0;

	  for (int ii=0;ii<N;ii++)
	    {
	      xcmd+=particleData(ii,id);
	    }
	  
	  positions(i,id)=xcmd;
	}
      
      storer.store(w,waveT);
      sP.accumulate(w,waveT);
      
    }
  
  auto & storedCorrelators = w.getStorageScalarCorrelators()[storerLabel];

  const auto & i = w.getTimeIndex().at(storerLabel);
  
  auto j =  wrapIndex(i  -1 ,recordSteps ) ;
  
  //std::cout << j << std::endl;

  
  for (int id=0;id<getDimensions();id++)
    {
      ASSERT_NEAR(positions(nStores-1,id) ,storedCorrelators(j*D + id) , 1e-4 );
    }
  
  
}



TEST(initial_configurations,read_directory)
{
  if (getDimensions() == 1)
  {
  auto configurations= readStatesFromDirectory("../tst/testInitalConfigurationVMC");

  int nWwalkers =  configurations.size() ;

  ASSERT_EQ(nWwalkers,1);

  std::vector<double> xConf
  {
    0.4334811868529716,0.5123232645232533,0.7660927469332877,-0.5032132098187615,0.17076099960880708,-0.15155238913361904,0.5368818932102931,0.1496488524427948,-0.46513468566258886,-0.410402300451225
  };
  
  const auto & xConfLoaded = configurations[0][0];

  for (size_t i=0;i<xConf.size();i++)
    {

      ASSERT_NEAR ( xConf[i] ,xConfLoaded(i) , 1e-6 ) << std::endl;
    }


  }



} 