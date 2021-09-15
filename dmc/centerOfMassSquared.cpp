#include "centerOfMassSquared.h"
#include "walkers.h"
#include "wavefunction/productWavefunction.h"
#include "geometry.h"

centerOfMassSquared::centerOfMassSquared(const json_t & j) : centerOfMassSquared::centerOfMassSquared(j["set"].get<int>() ) {}

real_t centerOfMassSquared::operator()(walker_t & w,wavefunction_t & psi)
{
  real_t cm=0;
  auto & positions = w.getStates()[setA];
  auto & geo=psi.getGeometry();
  diff= geo.differencesOneBody(positions,{0.,0.,0.});
  dis= norm(diff);
  cm=dis.sum()*dis.sum()/dis.size();
  return cm;
  
}
