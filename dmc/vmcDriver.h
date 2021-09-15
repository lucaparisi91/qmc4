#ifndef VMC_DRIVER_H
#define VMC_DRIVER_H
#include "traits.h"
#include <string>
#include <vector>
#include <deque>
#include "driver.h"


class walker;

class configurationsSaver
{
public:
  configurationsSaver(int nConfSave_=1);
  virtual void dump(const walker & w,int i);
  const auto & getFolder() const {return baseDir;}

  auto & nConfigurationsToSave() {return _nConfigurationsToSave;}
protected:
  
private:
  std::string baseDir;
  std::deque<json_t> jSonConfigurations;
  int _nConfigurationsToSave;
 
};


class vmcDriver : public driver
{
public:
	using walker_t = walker;

	vmcDriver(wavefunction_t * wave,real_t sigma_);
	
	void run(states_t & states,size_t nSamples);
	void step(); // perform a mc step
	void out(); // output summery of the blocks

	virtual void accumulate() ; // accumulate measurements

	auto & currentWalker() {return current_walker;}
	auto & oldWalker(){return old_walker;}
  auto & nConfigurationsToSave() {return confSave.nConfigurationsToSave(); }

  auto & getConfigurationsSaver() {return confSave;}
  
private:
  walker_t current_walker;
  walker_t old_walker;
  walker_t tmp_walker;
  metropolis metropolisObj;
  std::unique_ptr<mover> vmcMove;
  configurationsSaver confSave;
  
};



bool update(walker & w,productWavefunction & psi);





#endif
