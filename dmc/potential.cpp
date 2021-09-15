#include "potential.h"
#include "geometry.h"
#include "tableDistances.h"
#include <iostream>
#include "tools.h"
#include "walkers.h"


potential::potential(const geometry_t & geo_ ) : geo(&geo_)
{

};

potential1b::potential1b(const geometry_t & geo_ , int setA_) : potential::potential(geo_),_setA(setA_)
{
	
};

harmonicPotential::harmonicPotential(const geometry_t & geo,real_t freq , int setA ) : potential1b(geo,setA) ,omega(freq) {}

real_t harmonicPotential::operator()(const walker_t & w)
{
  const auto & dis= w.getTableDistances().distances(setA() );

  real_t sum=0;
  
  for(int i=0;i<dis.size();i++)
    {
      sum+=0.5 * omega * omega * dis(i) * dis(i);
    }  
  return sum;
};

sumPotentials::sumPotentials(std::vector<potential*> potentials_) :
  _potentials(potentials_)
{
  
}

real_t sumPotentials::operator()(const walker_t & state )
{
  real_t sum=0;
  for (int i=0;i<size();i++)
    {
      sum+=(*(_potentials[i]))(state);
    }
  return sum;
}

harmonicPotential::harmonicPotential(const json_t & j, const geometry_t & geo) : harmonicPotential(geo, j["omega"].get<real_t>(), j["set"].get<int>() ) {};

squareWellPotential2b::squareWellPotential2b(const geometry_t & geo,real_t V0_ , real_t R0_ , int setA_, int setB_ ) : R0(R0_),setA(setA_),setB(setB_),V0(V0_), potential(geo) {}

squareWellPotential2b::squareWellPotential2b(const json_t & j, const geometry_t & geo) : squareWellPotential2b(geo,j["V0"].get<real_t>(),j["R0"].get<real_t>(),j["sets"][0].get<int>(),j["sets"][1].get<int>()) {}


real_t squareWellPotential2b::operator()(const walker_t & w)
{
  const auto & dis= w.getTableDistances().distances(setA,setB);
  
  auto nPairsInt=std::count_if(dis.data(),dis.data() + dis.size(), [&] (const real_t  &d) {return d<= R0;});
  
  return -nPairsInt*V0;
  
};

dipolarPotential2b::dipolarPotential2b(const geometry_t & geo,real_t D_, int setA_, int setB_) : setA(setA_),setB(setB_),D(D_),k(0),potential(geo) {}


dipolarPotential2b::dipolarPotential2b(const json_t & j,const geometry_t & geo) : dipolarPotential2b::dipolarPotential2b(geo,j["D"].get<real_t>(),j["sets"][0].get<int>(), j["sets"][1].get<int>() ) {};

real_t dipolarPotential2b::operator()(const walker_t & w)
{
  double v=0;
  double d=0;

  const auto & norms= w.getTableDistances().distances(setA,setB);
  
  double lBox=getGeometry().getLBox(0);
  
  for (auto i=0;i<norms.size();i++) 
      {
	double x = norms(i);
	v+=D/(x*x*x);
	
	for (int n=1;n<=k;n++)
	  {
	    d=x + n*lBox;
	    v+=D/(d*d*d);
	    d= n*lBox-x;
	    v+=D/(d*d*d);
	    
	  }
      }
    
    return v;

}

real_t poschTellerPotential2b::operator()(const walker_t & w)
{
  
    real_t V=0;
    
    const auto & norms= w.getTableDistances().distances(setA,setB);
    
    for (auto i=0;i<norms.size();i++)
      {
	real_t tmp = cosh(norms(i)*alpha);
	V-=2 * alpha * alpha /( tmp*tmp);
      }
    return V;
    
}

poschTellerPotential2b::poschTellerPotential2b(const geometry_t & geo,real_t R0_, int setA_,int setB_) :
  potential(geo),
  R0(R0_),
  setA(setA_),
  setB(setB_),
  alpha(1/R0)
{
  
}

poschTellerPotential2b::poschTellerPotential2b(const json_t & j, const geometry_t & geo)
  : poschTellerPotential2b::poschTellerPotential2b(geo,j["R0"].get<real_t>(),j["sets"][0].get<int>(),j["sets"][1].get<int>())
{

}
