#include "traits.h"

class wavefunction;
class tableDistances;
class productWavefunction;
class walker;
class dmcWalker;
class energy;
class tableSlaters;


template<class T> class walkerContainer;

class estimatorCollection;
struct initializer
{
  static void registerDistances(tableDistances & tab,const wavefunction & wave);
  static void registerDistances(tableDistances & tab,const productWavefunction & wave);
  static void addDistances(walker & w,estimatorCollection & collection);
  
  static bool initialize(walker & w, const states_t & states ,  productWavefunction & psi);
  static bool initialize(dmcWalker & w, const states_t & states ,  productWavefunction & psi,energy & ob);

  static bool initialize(walkerContainer<dmcWalker> & ws, const std::vector<states_t> & states ,  productWavefunction & psi,energy & ob);

  static void registerSlaterOrbitals(tableSlaters & tab,const productWavefunction & psi);
  
  
};
