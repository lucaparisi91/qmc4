#include "orbitals.h"
#include "tools.h"

#if DIMENSIONS == 3
sinOrbital::sinOrbital(int n1,int n2,int n3,double lBox_) : lBox(lBox_) , ns{n1,n2,n3}
{
  
  k.resize(getDimensions() );
  k[0]=n1*2*M_PI/lBox;
  k[1]=n2*2*M_PI/lBox;
  k[2]=n3*2*M_PI/lBox;
  
  int isCos=1;
    
    if (abs(n1)>0 )
      {
	isCos=n1/std::abs(n1);
      }
    else
      {
	if (abs(n2)>0)
	  {
	    isCos=n2/std::abs(n2);
	  }
	else
	  {
	    if (abs(n3)>0)
	      {
		isCos=n3/std::abs(n3);
	      }
	  }
	  
      }
	
    if (isCos==1)
      {
	delta=M_PI/2.;
      }
    else
      {
	delta=0;
      }
   
  }

sinOrbital::sinOrbital(int n1,int n2,int n3,const json_t & j) : sinOrbital(n1,n2,n3,j["lBox"].get<real_t>()) {}

template<class orbital_t>
void orbitalSet<orbital_t>::storeEvaluate(const state_t & states,orbitalSet<orbital_t>::matrix_t & orbitalMatrix) const
{
  constexpr int D = getDimensions();
  const int N = getN(states);
  
  orbitalMatrix.resize(N,N);
  assert(orbitals.size() >= N);
  for (int j=0;j<N;j++)
    {
      auto & orbital = orbitals[j];
      for (int i=0;i<N;i++)
	{
	  orbitalMatrix(i,j)=orbital(    states(i,0) , states(i,1),states(i,2)   );	  
	}
    }
};


template<class orbital_t>
 orbitalSet<orbital_t>::orbitalSet(const json_t & j)
{
  int n=j["n"].get<int>();
  
  fillFermiSeaByEnergy(getOrbitals() , n, j );
  
}


planeWave::planeWave(int nx,int ny,int nz,real_t lBox,real_t teta)
{
  k.resize(getDimensions() );
  k[0]=(nx*2*M_PI + teta)/lBox ;
  k[1]=(ny*2*M_PI+teta)/lBox ;
  k[2]=(nz*2*M_PI + teta)/lBox ;
}



planeWave::planeWave(int nx,int ny,int nz,const json_t & jI) : planeWave(nx,ny,nz,jI["lBox"].get<real_t>(),    jI.find(std::string("teta") )!=jI.end() ? jI["teta"].get<real_t>() : 0  ) {}



template class orbitalSet<sinOrbital>;
template class orbitalSet<planeWave>;



#endif
