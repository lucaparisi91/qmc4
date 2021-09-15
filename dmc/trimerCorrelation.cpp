#include "trimerCorrelation.h"
#include "walkers.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include "tools.h"
#include "tools3B.h"


void trimerCorrelation::accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)
{
  assert(getN(w.getStates()[setC] )== 1);

  setNormalizationFactor(w,wavefunction,acc);
    auto & norms1 =w.getTableDistances().distances(setA,setC);
    auto & norms2 =w.getTableDistances().distances(setB,setC);

    auto & norms3 = w.getTableDistances().distances(setA,setB);
    real_t radius=0;

    auto N1 = norms1.size();
    auto N2 = norms2.size();
    auto N3 = norms2.size();


     
  for (int i=0;i<N1;i++ )
    for (int j=0;j<N2;j++ )
      {
	auto d3 =norms3(i*N2 + j);
	radius=std::sqrt( norms1(i)*norms1(i) + norms2(j)*norms2(j) + d3*d3 );
	if ( radius < acc.maxx() )
	  {
#if DIMENSIONS == 1
	      acc.accumulate(_normalizationFactor,radius);
#endif
#if DIMENSIONS == 3
	      
	    acc.accumulate(_normalizationFactor/(radius*radius),radius);
#endif
	  }
      }
  
  acc.weight()+=1;
}


trimerCorrelation::trimerCorrelation(int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_),_normalizationFactor(0)
{
    assert(setA!=setB);
    assert(setA!=setC);
}

trimerCorrelation::trimerCorrelation(const json_t & j) : trimerCorrelation(j["sets"][0].get<int>(),j["sets"][1].get<int>(),j["sets"][2].get<int>())
{
  
}



void trimerCorrelation::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const trimerCorrelation::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  auto NB = getN(w.getStates()[setB]);
  auto NC = getN(w.getStates()[setB]);
  
#if DIMENSIONS == 3
  _normalizationFactor=1/(dx*4*M_PI*NA*NB*NC);
#endif

#if DIMENSIONS == 1
  
  _normalizationFactor=std::pow(lBox,1)/(dx*NA*NB*NC*2);   
#endif
  
  
}



void trimerCorrelationUnDis::accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)
{

  setNormalizationFactor(w,wavefunction,acc);

  const auto & dis=w.getTableDistances().distances(setA,setA);
  real_t sum=0;

    const int N = getN(w.getStates()[setA]);

    
    LOOP3B( N,
	    auto rji = dis(ji);	
	    auto rkj = dis(kj);	
	    auto rki = dis(ki);
	    
	    auto radius= sqrt(rji*rji + rki*rki + rkj*rkj);

	    if ( radius < acc.maxx() )
	      {
#if DIMENSIONS == 1
		acc.accumulate(_normalizationFactor,radius);
#endif
#if DIMENSIONS == 3
	      
	    acc.accumulate(_normalizationFactor/(radius*radius),radius);
#endif
	      }
	    
	    
	    )
      
      acc.weight()+=1;
}


trimerCorrelationUnDis::trimerCorrelationUnDis(int setA_) : setA(setA_),_normalizationFactor(0)
{
}

trimerCorrelationUnDis::trimerCorrelationUnDis(const json_t & j) : trimerCorrelationUnDis(j["sets"][0].get<int>() )
{
  
}


void trimerCorrelationUnDis::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const trimerCorrelation::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  
#if DIMENSIONS == 3
  _normalizationFactor=1/(dx*4*M_PI*NA*NA*NA);
#endif

#if DIMENSIONS == 1
  _normalizationFactor=std::pow(lBox,1)/(dx*NA*NA*NA*2);   
#endif
  
  
}












void trimerCorrelationDis::accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)
{

  setNormalizationFactor(w,wavefunction,acc);

  real_t sum=0;

  const int NA = getN(w.getStates()[setA]);
  const int NB = getN(w.getStates()[setB]);
  const int NC = getN(w.getStates()[setC]);

  const auto & disAB = w.getTableDistances().distances(setA,setB);
  const auto & disBC = w.getTableDistances().distances(setB,setC);
  const auto & disAC = w.getTableDistances().distances(setA,setC);
  
    
  LOOP3B_DIS( NA, NB,NC , 
	      
	      auto rij = disAB(ij);	
	      auto rjk = disBC(jk);	
	      auto rik = disAC(ik);
	    
	      auto radius= sqrt(rij*rij + rik*rik + rjk*rjk);

	      if ( radius < acc.maxx() )
		{
#if DIMENSIONS == 1
		  acc.accumulate(_normalizationFactor,radius);
#endif
#if DIMENSIONS == 3
	      
		  acc.accumulate(_normalizationFactor/(radius*radius),radius);
#endif
		}
	    
	    
	      )
      
    acc.weight()+=1;
}


trimerCorrelationDis::trimerCorrelationDis(int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_) , _normalizationFactor(0)
{
  assert( (setA != setB) and (setB != setC) );
}

trimerCorrelationDis::trimerCorrelationDis(const json_t & j) : trimerCorrelationDis(j["sets"][0].get<int>() ,j["sets"][1].get<int>() ,j["sets"][2].get<int>() )
{
  
}


void trimerCorrelationDis::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const trimerCorrelationDis::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  auto  NB = getN(w.getStates()[setB]);
  auto  NC = getN(w.getStates()[setC]);
  
  
#if DIMENSIONS == 3
  _normalizationFactor=1/(dx*4*M_PI*NA*NB*NC);
#endif

#if DIMENSIONS == 1
  _normalizationFactor=std::pow(lBox,1)/(dx*NA*NB*NC*2);   
#endif
  
  
}
