#ifndef JASTROW_H
#define JASTROW_H

#include <vector>
#include <iostream>
#include "qmcExceptions.h"
#include <nlohmann/json.hpp>
#include "traits.h"
#include "parameters.h"

class optimizationParameter;
class mappedOptimizationParameter;


template<class T>
class jastrow
{
  /*
    CRT base class for the jastrow functions. Add default funcionalities to derived classes. 
    Derived class must posses at least a d0(x),d1(x)and d2(x) twhiche compute the zero,first and second derivative of the jastrow function
	*/
public:
  
  using gradient_t = std::vector<real_t>;
  
  using iterator_t = decltype( std::declval<gradient_t>().begin()  );

  
  
  void evaluateDerivatives(real_t x, real_t & d0_,real_t & d1_, real_t & d2_) {d0_=static_cast<const T*>(this)->d0(x);d1_=static_cast<const T*>(this)->d1(x);d2_=static_cast<const T*>(this)->d2(x);}
  
  void addGradientParameter(real_t x, const optimizationParameter & param, iterator_t begin, iterator_t end ) {throw missingImplementation("Jastrow does not seem to support any parameter gradient.");} ;

  void addGradientParameter(real_t x, const mappedOptimizationParameter & param, gradient_t & grad)
  {
    static_cast<T *>(this)->addGradientParameter(x,param,param.begin(grad),param.end(grad) );
  }
  

  int nParameters() {return 0;}
  
  std::string print(real_t minx,real_t maxx,size_t n)   
  {
    std::stringstream ss;
    assert(n>=1);
    
    real_t deltax=(maxx - minx)/(n);
    real_t d0,d1,d2;
    
    for(int i=0;i<n;i++)
      {
	real_t x = minx +  (i+0.5)*deltax;
	
	static_cast<T*>(this)->evaluateDerivatives(x,d0,d1,d2);
	
	ss << x << " " <<  d0 << " " << d1 << " "<< d2 << std::endl ;
	
      }

    
    return ss.str();
    
  }
  
private:

  

protected:
	jastrow(){}; // disalloes the instantation of base class. Only concrete jastrows shoulf derive from this class
  
};



#endif
