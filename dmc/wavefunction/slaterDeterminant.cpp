#include "slaterDeterminant.h"
#include "slaters.h"
#include "orbitals.h"
#include "walkers.h"

template<class orbitalSet_t>
slaterDeterminantWavefunction<orbitalSet_t>::slaterDeterminantWavefunction(const orbitalSet_t & orbitals_,const geometry_t & geo_,size_t setA_) : _orbitals( std::make_unique<orbitalSet_t>(orbitals_) ),setA(setA_),wavefunction::wavefunction(geo_) {}


template<class orbitalSet_t>
real_t  slaterDeterminantWavefunction<orbitalSet_t>::operator()(const walker_t & w)
{
  return w.getTableSlaters().logDeterminantReal(setA);
};

template<class orbitalSet_t>
void slaterDeterminantWavefunction<orbitalSet_t>::accumulateDerivatives(walker_t & w)
{
  
  const auto & inverseSlater = w.getTableSlaters().slaterMatrixRealInverse(setA);
  const auto & state = w.getStates()[setA];
  const auto & orbitals= _orbitals->getOrbitals();
  auto & lap=w.getLaplacianLog();
  
  auto & grad = w.getGradients()[setA];
  
  const auto & slater = w.getTableSlaters().slaterMatrixReal(setA);
  real_t tmp2=0;
  
  const int N = getN(state);
  
  constexpr int D = getDimensions();
  std::array<real_t,3> dx;
  assert( N <= orbitals.size());

  tmpGrad.resize(N,D);
  tmpGrad.setConstant(0);
  for(int j=0;j<N;j++)
      {
	const auto & orbital = orbitals[j];
	
	for(size_t i=0;i<N;i++)
	  {
	    
	    orbital.evaluateDerivatives( state(i,0),state(i,1),state(i,2) , dx[0],dx[1],dx[2],tmp2);
	    for(size_t id=0;id<D;id++)
	      {
		tmpGrad(i,id)+=dx[id]*inverseSlater(j,i);
	      }
	    
	    lap+=tmp2*inverseSlater(j,i);
	    
	  }
	    

	    
		    	  
      }
  lap-= (tmpGrad*tmpGrad).sum();
  grad+=tmpGrad;
  w.getLogWave()+=(*this)(w);
}

template<class orbitalSet_t>
slaterDeterminantWavefunction<orbitalSet_t>::slaterDeterminantWavefunction(const json_t & j,const geometry_t & geo) : slaterDeterminantWavefunction( orbitalSet_t(j["orbitalSet"] ),geo,j["set"].get<int>()) {}


template<class orbitalSet_t>
slaterDeterminantWavefunctionComplex<orbitalSet_t>::slaterDeterminantWavefunctionComplex(const orbitalSet_t & orbitals_,const geometry_t & geo_,size_t setA_) : _orbitals( std::make_unique<orbitalSet_t>(orbitals_) ),setA(setA_),wavefunction::wavefunction(geo_) {}



template<class orbitalSet_t>
real_t  slaterDeterminantWavefunctionComplex<orbitalSet_t>::operator()(const walker_t & w)
{
  return w.getTableSlaters().logDeterminantComplex(setA).real();
};



template<class orbitalSet_t>
void slaterDeterminantWavefunctionComplex<orbitalSet_t>::accumulateDerivatives(walker_t & w)
{
  using value_t = std::complex<real_t> ;
  const auto & inverseSlater = w.getTableSlaters().slaterMatrixComplexInverse(setA);
  const auto & state = w.getStates()[setA];
  const auto & orbitals= _orbitals->getOrbitals();
  value_t lap=0;
  
  auto & grad = w.getGradients()[setA];
  auto & gradPhase = w.getPhaseGradients()[setA];
  
  const auto & slater = w.getTableSlaters().slaterMatrixComplex(setA);
  std::complex<real_t> tmp2=0;
  
  const int N = getN(state);
  
  constexpr int D = getDimensions();
  std::array<value_t,3> dx;
  
  assert( N <= orbitals.size());

  
  tmpGrad.resize(N,D);
  tmpGrad.setConstant(0);
  for(int j=0;j<N;j++)
      {
	const auto & orbital = orbitals[j];
	
	for(size_t i=0;i<N;i++)
	  {
	    
	    orbital.evaluateDerivatives( state(i,0),state(i,1),state(i,2) , dx[0],dx[1],dx[2],tmp2);
	    for(size_t id=0;id<D;id++)
	      {
		
		tmpGrad(i,id)+=dx[id]*inverseSlater(j,i);
	      }
	    
	    lap+=(tmp2*inverseSlater(j,i));
	    
	  }
	
      }

  
 
  for(int i=0;i<grad.size();i++)
    {
      lap-=(tmpGrad(i)*tmpGrad(i));
      grad(i)+=tmpGrad(i).real();
      gradPhase(i)+=tmpGrad(i).imag();
    }
  
  w.getLogWave()+=(*this)(w);
  w.getLaplacianLog()+=lap.real();

}

template<class orbitalSet_t>
slaterDeterminantWavefunctionComplex<orbitalSet_t>::slaterDeterminantWavefunctionComplex(const json_t & j,const geometry_t & geo) : slaterDeterminantWavefunctionComplex( orbitalSet_t(j["orbitalSet"] ),geo,j["set"].get<int>()) {}



#if DIMENSIONS == 3

template class slaterDeterminantWavefunction<orbitalSet<sinOrbital> > ;
template class slaterDeterminantWavefunctionComplex<orbitalSet<planeWave> > ;

#endif
