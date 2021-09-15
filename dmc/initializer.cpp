#include "initializer.h"
#include "wavefunction/productWavefunction.h"
#include "tableDistances.h"
#include "qmcExceptions.h"
#include "walkers.h"
#include "wavefunction/wavefunction.h"
#include "estimatorCollection.h"
#include "estimators.h"

void initializer::addDistances(walker & w,estimatorCollection & collection)
{
  for (auto & est : collection)
    {
      w.getTableDistances().add(est->sets());
    }
  
  for (const auto & storeO : collection.storers())
    {
      w.getTableDistances().add(storeO->sets());
    }
  
}
 

void initializer::registerDistances(tableDistances & tab,const wavefunction & wave)
{
	const auto & sets = wave.sets();
	tab.setGeometry(wave.getGeometry());
	tab.add(sets);

};

void initializer::registerDistances(tableDistances & tab,const productWavefunction & waves)
{
	
	for (const auto & wave : waves.waves() )
		registerDistances(tab,*wave);
};

void initializer::registerSlaterOrbitals(tableSlaters & tab,const productWavefunction & psi)
{
  for (const auto & wave : psi.waves() )
    {
      const auto  sets = wave->sets();
      const auto  orbitals = wave->orbitals();

      if (orbitals.size() == 1 and sets.size() == 1 )
	{
	  tab.add(sets[0],orbitals[0]);
	}
    }
}

bool initializer::initialize(walker & w, const states_t & states ,  productWavefunction & psi)
{
	w.getStates()=states;
	
	registerDistances(w.getTableDistances(),psi);
	registerSlaterOrbitals(w.getTableSlaters(),psi);
	w.getTableDistances().update(states);
	w.getTableSlaters().update(states);
	if ( psi.satisfyConstraints(w) )
	  {
	    psi.evaluateDerivatives( w);
	    return true;
	  }
	else
	  {
	    return false;
	  }
	
	
}

bool initializer::initialize(dmcWalker & w, const states_t & states ,  productWavefunction & psi,energy & ob)
{
  w.getStates()=states;
  registerDistances(w.getTableDistances(),psi);
  registerSlaterOrbitals(w.getTableSlaters(),psi);
  return updateForceGradientEnergy(w,psi,ob);
  
}

bool initializer::initialize(walkerContainer<dmcWalker> & ws, const std::vector<states_t> & states ,  productWavefunction & psi,energy & ob)
{
  
  ws.resize(states.size());
  
  for (int i=0;i<states.size();i++)
    {
      bool status = initialize(ws[i],states[i],psi,ob);
      if (!status)
	{
	  return false;
	}
    }

  return true;
}


  
