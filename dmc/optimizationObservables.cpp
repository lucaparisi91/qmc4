#include "optimizationObservables.h"
#include "energy.h"
#include "wavefunction/productWavefunction.h"
#include "parameters.h"



optimizationObservables::optimizationObservables(int nP,std::vector<std::vector<mappedOptimizationParameter> > parameters_ , energy & energyOb_) :
  realVectorObservable::realVectorObservable(2*nP + nP*(nP+1)/2),
  energyOb(&energyOb_),
  parameters(parameters_)
{
  parameterGradient.resize(nP);
  
}

void optimizationObservables::accumulate(walker_t & w, wavefunction_t & wave,accumulator_t & acc)
{
  std::fill(parameterGradient.begin(),parameterGradient.end() , 0 );
  real_t e=(*energyOb)(w,wave);
  
  wave.addGradientParameter(w,parameters,parameterGradient);
  
  int Np = parameterGradient.size();
  real_t tmp=0;
  
  size_t p=0;
  
  for(int i=0;i<Np;i++)
    {
      acc.accumulate(i,parameterGradient[i] );
      p++;
    };
  
  for(int i=0;i<Np;i++)
    {
      tmp=parameterGradient[i]*e;
      acc.accumulate(p,tmp);
      p++;
    }
  
  for(int i=0;i<Np;i++)
    {
      for(int j=0;j<=i;j++)
	{
	  tmp=parameterGradient[i]*parameterGradient[j];
	  acc.accumulate(p,tmp);
	  p++;
	}
    }
  
  acc.accumulate(p,e);
  
}
