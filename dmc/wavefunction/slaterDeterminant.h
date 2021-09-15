#include "wavefunction.h"

template<class orbitalSet_t>
class slaterDeterminantWavefunction : public wavefunction
{
public:
  using grad_t = ::difference_t;
  
  slaterDeterminantWavefunction(const orbitalSet_t &  orbital, const geometry_t & geo,size_t setA);
  slaterDeterminantWavefunction(const json_t & j,const geometry_t & geo);
  
  
  virtual void accumulateDerivatives(walker_t & w);

  
  virtual real_t operator()(const walker_t & w);

  virtual std::vector<int> sets() const override {return {setA};}
  
  std::vector<orbitalSetBase*> orbitals() const override {return {_orbitals.get()};}

  static std::string name() {return "slater/" + orbitalSet_t::name();}
  
private:
  grad_t tmpGrad;
  std::unique_ptr<orbitalSet_t> _orbitals;
  int setA;
};


template<class orbitalSet_t>
class slaterDeterminantWavefunctionComplex : public wavefunction
{
public:
  using grad_t = ::difference_complex_t;
  
  slaterDeterminantWavefunctionComplex(const orbitalSet_t &  orbital, const geometry_t & geo,size_t setA);
  slaterDeterminantWavefunctionComplex(const json_t & j,const geometry_t & geo);
  
  
  virtual void accumulateDerivatives(walker_t & w);

  virtual bool isComplex() const override {return true;}
  
  virtual real_t operator()(const walker_t & w);
  
  virtual std::vector<int> sets() const override {return {setA};}
  
  std::vector<orbitalSetBase*> orbitals() const override {return {_orbitals.get()};}
  
  static std::string name() {return "slaterComplex/" + orbitalSet_t::name();}
  
private:
  grad_t tmpGrad;
  std::unique_ptr<orbitalSet_t> _orbitals;
  int setA;
};


