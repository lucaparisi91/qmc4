#ifndef ESTIMATORS_COLLECTION_H
#define ESTIMATORS_COLLECTION_H


#include <vector>
#include <memory>
#include "traits.h"
#include "correlationEstimator.h"

class walker;
class estimatorBase;
class storer;


class estimatorCollection
{
public:
  using walker_t =  walker;
  using wavefunction_t = ::wavefunction_t;
  
  estimatorCollection(){};
  auto begin() {return _estimators.begin();}
  const auto cbegin() const {return _estimators.begin();}

  auto end() {return _estimators.end();}
  const auto cend() const {return _estimators.end();}

  void push_back(storer* storer_){_storers.emplace_back(storer_); }

  
  void push_back(estimatorBase* est){_estimators.emplace_back(est);}
  

  void accumulate(walker_t & w,wavefunction_t & psi);
  
  auto & operator[](size_t i) {return _estimators[i];}
  virtual void clear();
  virtual void dump();
  virtual void accumulateMPI(int root=0);
  virtual void reserve(walker_t & w);

  const auto & storers() const {return _storers;}
  
  
private:
  std::vector<std::unique_ptr<estimatorBase> > _estimators;
  std::vector<std::unique_ptr<storer> > _storers;
  
};

#endif
