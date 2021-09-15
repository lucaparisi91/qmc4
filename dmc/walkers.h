#ifndef WALKERS_H
#define WALKERS_H

#include "traits.h"
#include <memory>
#include "tableDistances.h"
#include "qmcExceptions.h"
#include "slaters.h"
#include "ptools.h"
#include <fstream>

/*
A walker contains all the informiation
on the current configurations. 
Also maintains cached data for use, as particle distances , wavefunction values and
so on. Walkers own data memory
*/


struct walker
{
public:
  using grads_t = ::states_t;
  walker(){};

  walker(const walker & w2 ) ;

  void operator=(const walker & w2 );
  const  auto & getStates() const {return _states;}
  const auto & getTableDistances() const {return _tab;}
  const auto & getLogWave() const {return _waveValue;}
  const auto & getGradients() const {return _gradients;}
  const auto & getLaplacianLog() const {return _lapLog;}
  const auto & getTableSlaters()  const {return _slaters;}


  auto & getStates()  {return _states;}
  auto & getTableDistances()  {return _tab;}
  auto & getTableSlaters()  {return _slaters;}
  
  auto & getLogWave() {return _waveValue;}
  auto & getGradients()  {return _gradients;}
  auto & getPhaseGradients()  {return _phaseGradients;}
  auto & getLaplacianLog() {return _lapLog;}

  auto & getStorageScalarCorrelators() {return storageScalarObservables;}
  const auto & getStorageScalarCorrelators() const {return storageScalarObservables;}
  
  const auto & getTimeIndex() const {return correlatorCurrentTimeIndex;}
  auto & getTimeIndex()  {return correlatorCurrentTimeIndex;}
  
  
  virtual const real_t & getEnergy() const {throw missingImplementation("Energy not accessible from the walker"); return _waveValue;};
  virtual real_t & getEnergy()  {throw missingImplementation("Energy not accessible from the walker");return _waveValue;};

  auto &  getMPIDatatype() {return dtype;}
  
  virtual void updateMPIDataType()  {throw missingImplementation("Creation of MPI walker data type");};

  auto & getFillingStatus() {return isFilling;}
  const auto & getFillingStatus() const {return isFilling;}
  
private:
  states_t _states; // a vector of particle data
  tableDistances _tab;
  real_t _waveValue; // value of the wavefunction
  grads_t _gradients; // contains the gradient of the wavefunction. Just a temporary
  real_t _lapLog; // contains the laplacian of the logarithm of the wavefunction
  tableSlaters _slaters; // contains the matrix of slater determinants
  MPI_Datatype dtype;
  grads_t _phaseGradients; // contains the gradient of the wavefunction. Just a temporary
  std::map<std::string,Eigen::VectorXd> storageScalarObservables;
  
  

  std::map<std::string,int> correlatorCurrentTimeIndex;
  std::map<std::string,bool> isFilling;
};

struct dmcWalker : public walker
{
  dmcWalker();
  dmcWalker(const dmcWalker & w);
  void operator=(const dmcWalker & w);
  
  
  virtual real_t & getEnergy() override {return _e;}
  virtual const real_t & getEnergy() const override {return _e;}
  virtual void updateMPIDataType() ;
  
private:

  void createMPIDataType() ;
  
  real_t _e=1E+20; // stores the energy of the current configuration
 
};




class energy;


void update(dmcWalker & w, productWavefunction & wave);
bool updateForceGradientLaplacian(walker & w,productWavefunction & psi);
bool updateForceGradientEnergy(dmcWalker & w,productWavefunction & psi, energy & energyOb);

template<class T>
class walkerContainer
{
public:
  using value_type=T;
  /*A wrapper around std::vector to contain walkers. 
    Calling resize does not make old allocated data invalid.
*/
  
  walkerContainer() : walkers(),_size(0),baseDir("configurations")
  {
    
  };
  //walkerContainer( std::vector<T> vec ) : walkers(vec),_size(vec.size()) {}
  
  auto & operator[](size_t i) {return *(walkers[i]);}
  const auto & operator[](size_t i) const {return *(walkers[i]);}
  
  void push_back(const T &  w);


 
  void resize(size_t size2);
  void resize(size_t size2, const T & w);
  
  void reserve(size_t size2,const T & w);
  size_t size() const {return _size;}
  size_t capacity() const {return walkers.size();}
  
  auto  begin()  {return walkers.begin();}
  auto end() {return walkers.begin() + _size;}

  auto  begin() const  {return walkers.begin();}
  auto end() const {return walkers.begin() + _size;}
  
  auto  cbegin() const {return walkers.cbegin();}
  auto cend() const {return walkers.cbegin() + _size;}

  auto data() {return walkers.data();}

  json_t toJson();
  
  void dump(int i);
  
  private:
  std::vector<std::unique_ptr<T> > walkers;
  size_t _size;
  std::string baseDir;
  bool saveOnlyLastConfiguration;
  
};


#endif
