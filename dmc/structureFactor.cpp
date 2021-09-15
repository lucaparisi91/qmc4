#include "structureFactor.h"
#include <complex>
#include "walkers.h"
#include "tools.h"

#if DIMENSIONS == 1

double getDeltaQ(double l_box,double qMax,int bins)
{
  
  double deltaQ;
  deltaQ=int(qMax/(2*M_PI*bins/l_box))*2*M_PI/l_box;
  
  if (deltaQ==0)
    {
      deltaQ=2*M_PI/l_box;
    }
  
  return deltaQ;
  
}

auto build_q_vector(int bins,double l_box,double qMax)
{
  structureFactor::vec_t qs(bins);
  double deltaQ;
  
  int i;
   
  deltaQ=getDeltaQ(l_box,qMax,bins);
   
  qs(0)=2*M_PI/l_box;
  
  for(i=1;i<bins;i++)
    {
      qs(i)=qs(0)+i*deltaQ;
    }
  return qs;
  
}

structureFactor::structureFactor(int setA,int setB,const structureFactor::vec_t & qs,int channel) : realVectorObservable::realVectorObservable(qs.size()),_setA(setA),_setB(setB),_qs(qs),_channel(channel)
{
  
}

structureFactor::structureFactor(const json_t & j) :
  structureFactor(
		  j["sets"][0].get<int>(),
		  j["sets"][1].get<int>(),
		  build_q_vector(j["bins"].get<int>(),j["lBox"].get<real_t>(),j["maxq"].get<real_t>() ),
		  j["channel"].get<int>()
		  )
{}

void structureFactor::accumulate(structureFactor::walker_t & w,structureFactor::wavefunction_t & wave,structureFactor::accumulator_t & acc)
  {
  
    auto & p1= w.getStates()[_setA];
    auto & p2= w.getStates()[_setB];
    
    std::complex<double> sF1(0,0);
    std::complex<double> sF2(0,0);
    
    for(unsigned int j=0;j<_qs.size();j++)
      {
	sF1=(0,0);
	for(int i=0;i<getN(p1);i++)
	  {
	    sF1+=std::exp(std::complex<double>( 0,  _qs[j]*p1(i,0) ) );
	  }
	if (_setA != _setB)
	  {
	    sF2=(0,0);
	    for(int i=0;i<getN(p2);i++)
	      {
		sF2+=std::exp(std::complex<double>( 0,  _qs[j]*p2(i,0) ) );
	      }
	    acc.accumulate(std::norm(sF1 + std::complex<double>(_channel,0) * sF2)/(getN(p1) + getN(p2) ) ,j);

	    
	  }
	else
	  {
	    acc.accumulate(std::norm(sF1)/(getN(p1) ) ,j);

	  }
	
      }
  }


#endif
