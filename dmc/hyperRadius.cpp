#include "hyperRadius.h"
#include "walkers.h"
#include "tools.h"
trimerhyperRadius::trimerhyperRadius(int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_)
{
  assert(setA!=setB);
  assert(setA!=setC);
  
}

trimerhyperRadius::trimerhyperRadius(const json_t & j) : trimerhyperRadius(j["sets"][0].get<int>(),j["sets"][1].get<int>(),j["sets"][2].get<int>())
{
  
  
}

real_t trimerhyperRadius::operator()(walker_t & w,wavefunction_t & psi)
{
  auto & norms1 =w.getTableDistances().distances(setA,setC);
  auto & norms2 =w.getTableDistances().distances(setB,setC);

  
  assert(  getN(w.getStates()[setC]) == 1);
  auto & norms3 = w.getTableDistances().distances(setA,setB);
  real_t radius=0;

  auto N1 = norms1.size();
  auto N2 = norms2.size();
  auto N3 = norms2.size();
  int k=0;
  
  for (int i=0;i<N1;i++ )
    for (int j=0;j<N2;j++ )
      {
	auto d3 =norms3(i*N2 + j);
	radius+=std::sqrt( norms1(i)*norms1(i) + norms2(j)*norms2(j) + d3*d3 );
      }

  
  radius=radius/(N1*N2);
  
  return radius;
  
}
