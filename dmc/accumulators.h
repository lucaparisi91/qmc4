#ifndef ACCUMULATORS_H
#define ACCUMULATORS_H

#include "traits.h"
#include "ptools.h"
#include <iostream>

template<class T>
class accumulator
{
public:
  void writeAverage(std::ostream & o) const
  {
    o << static_cast<const T*>(this)->average();
  }
};


template<class T>
class scalarAccumulator : public accumulator<scalarAccumulator<T> >
{
public:
  using value_t = T;
  scalarAccumulator() : sum(0.),n(0) {}

  void operator+=(value_t e){sum+=e;n+=1;};
  
  value_t average() const {return sum/n;}
  
  void clear(){sum=0.;n=0.;}
  
  void accumulateMPI(int root) // sum partial sums on all processors into the root processor
  {  
    sum=pTools::sum(sum,root);
    n=pTools::sum(n,root);
  }

  const auto & getWeight() {return n;}

 

  
private:
	value_t sum;
	size_t n;
};

template<class T>
class vectorAccumulator : public accumulator<vectorAccumulator<T> >
{
public:
  using vec_t = Eigen::Array<real_t,Eigen::Dynamic,1>;
  using value_t = T ;
  vectorAccumulator(){};
  vectorAccumulator(size_t size) : _sums(size),ns(size){
    _sums.setConstant(0);
    ns.setConstant(0);
  }
  
  void resize(size_t size, T val = 0)
  {
    _sums.resize(size);
    ns.resize(size);
    _sums.setConstant(val);
    ns.setConstant(val);
  }
    
  void accumulate(value_t a, size_t i,real_t weight=1) {_sums[i]+=a;ns[i]+=weight;}
  
  vec_t average() const {return _sums/ns;}

  const vec_t & sums() const {return _sums;}

  const vec_t & weights() const {return ns;}

  size_t size(){return  _sums.size();}

  void clear()
  {
    ns.setConstant(0);
    _sums.setConstant(0);
  }

  void accumulateMPI(int root)
  {
    pTools::sum(_sums,root);
    pTools::sum(ns,root);
  }
private:
  vec_t _sums;
  vec_t ns;
};


template<class T>
class histogramAccumulator : public accumulator<histogramAccumulator<T> >
{
public:
  using vec_t = Eigen::Matrix<T,Eigen::Dynamic,1> ;
  using value_t = T ;
  histogramAccumulator(){};

    void resize(size_t size,real_t min_,real_t max_)
  {
    _sums.resize(size);
    _minx=min_;
    _maxx=max_;
    _sums.setConstant(0);
    deltax=(_maxx-_minx)/size;
    deltaxInverse=1./deltax;
  }


  
  histogramAccumulator(size_t size,real_t min_, real_t max_) : _weight(0) {resize(size,min_,max_);}
  
  void accumulate(value_t a,value_t x) {int i=int((x-_minx)*deltaxInverse); _sums[i]+=a;}
  
  vec_t average() const {return sums()/weight();}
  
  const vec_t & sums() const {return _sums;}
   vec_t & sums()  {return _sums;}

  size_t size() const {return  _sums.size();}

  real_t maxx() const {return _maxx;}

  real_t & minx() {return _minx;}
  
  
  auto & weight() {return _weight;}
  const auto & weight() const {return _weight;}

  real_t stepSize() const {return deltax;}

  void clear(){_sums.setConstant(0.);_weight=0.;}
  
  void accumulateMPI(int root) // sum partial sums on all processors into the root processor
  {
    
    pTools::sum(_sums,root);
    _weight=pTools::sum(_weight,root);

  }

  auto x(int i) const {return _minx + (i+0.5)*deltax;}
  
  void writeAverage(std::ostream & o) const
  {
    auto av = average();

    for (int i=0;i<av.size();i++)
      {
	o << x(i) << " " << av(i) << std::endl;
      }
  }
  
private:
  vec_t _sums;
  T _weight;
  real_t _maxx;
  real_t _minx;
  real_t deltax;
  real_t deltaxInverse;
};



using realScalarAccumulator_t = scalarAccumulator<real_t>;
using realVectorAccumulator_t = vectorAccumulator<real_t>;
using realHistogramAccumulator_t = histogramAccumulator<real_t>;
#endif
