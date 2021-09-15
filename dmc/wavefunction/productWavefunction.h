#ifndef PRODUCT_WAVEFUNCTION_H
#define PRODUCT_WAVEFUNCTION_H

#include "traits.h"
#include <initializer_list>


class wavefunction;
class walker;
class mappedOptimizationParameter;


class productWavefunction
{
	/*
	Represent the sum of the logharitms of single wavefunctions. Does not own
	the wavefunctions.
	*/
public:
  
  using walker_t = walker;
  using gradientParameter_t=std::vector<real_t>;
  
  
  productWavefunction(){};
  

  productWavefunction(std::vector<wavefunction *> waves) ;

  real_t operator()(const walker_t & states); // evaluates \sum_i psi_i

  void evaluateDerivatives(walker_t & states); // computes the gradient and computes lap: log(\prod_i psi_i) and laplacian \nabla \psi = \prod_i exp(psi_i), lap_force : sum(grad**2) 

  void add(wavefunction * wave) ;


  size_t size () const { return _logWaves.size();}

  const auto & waves() const {return _logWaves;}

  const wavefunction & operator()(size_t i) const {return *(_logWaves[i]);}
  wavefunction & operator()(size_t i)  {return *(_logWaves[i]);}


  const geometry_t & getGeometry() const ;

  bool isComplex() {return isAtLeatOneComponentComplex;}

  bool satisfyConstraints(const walker_t & state);
  
  virtual void addGradientParameter(walker_t & w, const std::vector<std::vector< mappedOptimizationParameter> > & parameters,gradientParameter_t & grad  ) ;
  
  

  
private:
  std::vector<wavefunction*> _logWaves;
  bool isAtLeatOneComponentComplex;
};

#endif
