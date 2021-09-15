#include "wavefunction.h"
#include "tableDistances.h"
#include "parameters.h"

wavefunction::wavefunction(const geometry_t & geo_ ) : geo(&geo_)
{
  
};


bool wavefunction::satisfyConstraints(const walker_t & state)
{
  for (auto & storedConstraint : constraints)
    {
      if ( not (*storedConstraint)(state) )
	{
	  return false;
	}
    }
  return true;
};

void wavefunction::addConstraint( std::unique_ptr<constraint> p)
{
  constraints.push_back( std::move(p) );
}

void wavefunction::addGradientParameter(wavefunction::walker_t & w, const std::vector<mappedOptimizationParameter> & parameters, wavefunction::gradientParameter_t & grad )
  {
    for (auto & param : parameters )
      {
	addGradientParameter(w,param,grad);
      }
    
  }
  

void wavefunction::addGradientParameter(wavefunction::walker_t & w, const mappedOptimizationParameter & parameter, wavefunction::gradientParameter_t & grad )
  {
    auto it1 = parameter.begin(grad);
    auto it2 = parameter.end(grad);
    
    addGradientParameter(w,parameter,it1,it2);
    
  }
