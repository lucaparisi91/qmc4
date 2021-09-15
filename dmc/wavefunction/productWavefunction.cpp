#include "productWavefunction.h"
#include "tools.h"
#include "wavefunction.h"
#include "walkers.h"

void productWavefunction::evaluateDerivatives(productWavefunction::walker_t & w)
{
  auto & waveValue = w.getLogWave();
  auto & lap = w.getLaplacianLog();
  auto & grads = w.getGradients();
  auto & states = w.getStates();
  auto & gradsPhase = w.getPhaseGradients();
  
  lap=0;waveValue=0;
  grads.resize(states.size());
  
  if (isComplex() )
    {
      gradsPhase.resize(states.size());

      for (int i=0;i< gradsPhase.size();i++)
	{
	  gradsPhase[i].resize(getN(states[i]),getDimensions());
	  gradsPhase[i].setConstant(0.);
	}
  
    }
  for (int i=0;i< grads.size();i++)
    {
      grads[i].resize(getN(states[i]),getDimensions());
      grads[i].setConstant(0.);
    }
  
  for (int i=0;i< size();i++)
    {
      _logWaves[i]->accumulateDerivatives(w);
    }

}



real_t productWavefunction::operator()(const productWavefunction::walker_t &states)
{
  real_t waveValue=0;
  for (const auto & wave: _logWaves) 
    waveValue+=(*wave)(states);
  return waveValue;
}


const geometry_t & productWavefunction::getGeometry() const {return (_logWaves[0])->getGeometry();}

void productWavefunction::add(wavefunction * wave)
{
  _logWaves.push_back(wave);
  isAtLeatOneComponentComplex=isAtLeatOneComponentComplex or wave->isComplex() ;
  //std::cout << wave->isComplex() << std::endl;
}

productWavefunction::productWavefunction(std::vector<wavefunction *> waves) : productWavefunction()
{
  for ( auto & wave : waves )
    {
      add(wave);
    }
}


bool productWavefunction::satisfyConstraints(const walker_t & state)
{
  
  for (auto & wave : waves() )
    {
      if (not wave->satisfyConstraints(state) )
	{
	  return false;
	}
    }
  
  return true;
}


void productWavefunction::addGradientParameter(walker_t & w, const std::vector<std::vector< mappedOptimizationParameter> > & parameters,gradientParameter_t & grad  )
{
  for (int i=0;i<_logWaves.size();i++)
    {
      _logWaves[i]->addGradientParameter(w,parameters[i],grad);
    }
}
