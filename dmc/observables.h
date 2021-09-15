#ifndef OBSERVABLES_H
#define OBSERVABLES_H


#include "accumulators.h"
class productWavefunction;
class walker;

template<class T>
class observable
{
public:
	observable() {}
        using accumulator_t = T;
	using value_t = typename T::value_t ;
	using wavefunction_t = productWavefunction;
	using walker_t = walker;

        virtual void accumulate(walker_t & walker, wavefunction_t & wavefunction,  accumulator_t & acc)=0;

};
      

class realScalarObservable : public observable<realScalarAccumulator_t> 
{
public:
  virtual void accumulate(walker_t & w, wavefunction_t & wavefunction,  accumulator_t & acc) override
  {
    acc+=(*this)(w,wavefunction);
  }
  
  virtual real_t operator()(walker_t & walker,wavefunction_t & wave)=0;
  virtual std::vector<int> sets() const  {return {}; }
  
  static  bool isScalar() {return true;}

  static std::string kind() {return "scalar";}
  
};

class realHistogramObservable : public observable<realHistogramAccumulator_t>
{
public:
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)=0;
  virtual std::vector<int> sets() const  {return {}; }
  
  static bool isScalar() {return false;}
  static std::string kind() {return "histogram";}
  
};

class realVectorObservable : public observable<realVectorAccumulator_t>
{
public:
  using vec_t = Eigen::VectorXd;

  realVectorObservable(size_t size) : _size(size)
  {
    defaultx.resize(_size);
    for (int i=0;i<size;i++)
      {
	defaultx[i]=i;
      }
    
  }
  
  virtual void accumulate(walker_t & w,wavefunction_t & wavefunction,accumulator_t & acc)=0;
  virtual std::vector<int> sets() const  {return {}; }
  static bool isScalar() {return false;}
  virtual const vec_t & x() {return defaultx;}
  
  static std::string kind() {return "vector";}
  
private:
  real_t _size;
  vec_t defaultx;
  
  
};

#endif
