#include "pairCorrelation.h"
#include "tableDistances.h"
#include "walkers.h"
#include "tools.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"
#include <iostream>

pairCorrelation::pairCorrelation(int setA_,int setB_) : setA(setA_),setB(setB_),_normalizationFactor(0)
{}

pairCorrelation::pairCorrelation(const json_t & j) : pairCorrelation(j["sets"][0].get<int>(),j["sets"][1].get<int>())
{
  
}

void pairCorrelation::accumulate(walker_t & w,wavefunction_t & wave,pairCorrelation::accumulator_t & acc )
{
  
  auto & norms =w.getTableDistances().distances(setA,setB);
  
  setNormalizationFactor(w,wave,acc);
  for (int i=0;i<norms.rows();i++)
    {
      if (norms(i) < acc.maxx() )
  	{
#if DIMENSIONS == 3
  	  acc.accumulate(_normalizationFactor/(norms(i)*norms(i)),norms(i));
#endif
	  
#if DIMENSIONS == 1
  	  acc.accumulate(_normalizationFactor,norms(i));
#endif
	  
  	}
    }
  
  acc.weight()+=1;
}

void pairCorrelation::setNormalizationFactor(const walker_t & w , const wavefunction_t & psi ,const pairCorrelation::accumulator_t & acc) 
{
  auto lBox = psi.getGeometry().getLBox(0);
  auto dx = acc.stepSize();
  auto  NA = getN(w.getStates()[setA]);
  auto NB = getN(w.getStates()[setB]);
#if DIMENSIONS == 3
  _normalizationFactor=std::pow(lBox,3)/(dx*4*M_PI*NA*NB);
#endif

#if DIMENSIONS == 1
  
  _normalizationFactor=std::pow(lBox,1)/(dx*NA*NB*2);   
#endif
  
  _normalizationFactor = (setA == setB) ? 2*_normalizationFactor : _normalizationFactor;
  

}
