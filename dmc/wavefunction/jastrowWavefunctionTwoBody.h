#include "wavefunction/wavefunction.h"


template<class jastrow_t>
class jastrowTwoBodyWavefunctionDistinguishable : public wavefunction
{
public:
  
  jastrowTwoBodyWavefunctionDistinguishable(jastrow_t J_,const geometry_t  &geo_, int setA_=0,int setB_=0) : setA(setA_), setB(setB_),J(J_),wavefunction::wavefunction(geo_)
  {
    if (setA == setB) throw invalidInput("setA == setB in distinguishable two body jastrow");
  }
  
  jastrowTwoBodyWavefunctionDistinguishable(const json_t & j,const geometry_t & geo ) : jastrowTwoBodyWavefunctionDistinguishable( jastrow_t(j["jastrow"]), geo,j["sets"][0].get<int>() ,j["sets"][1].get<int>() ) {}
  
 virtual real_t operator()(const walker_t & walker) 
	{		
	  auto & dis = walker.getTableDistances().distances(setA,setB);
	  
	  int N=size(dis);
	  real_t sum=0;

	  for(int i=0;i<N;i++)
	    {
	      sum+=J.d0(dis(i));
	    }
	  return sum;
	};

  virtual std::vector<int> sets() const {return {setA,setB} ;}
  
  
  virtual void accumulateDerivatives( walker_t & walker ) override
  {
    auto & stateA = walker.getStates()[setA];
    auto & stateB = walker.getStates()[setB];

    auto & gradientA = walker.getGradients()[setA];
    auto & gradientB = walker.getGradients()[setB];
    
    auto & laplacian = walker.getLaplacianLog();
    auto & distances = walker.getTableDistances().distances(setA,setB);
    auto & differences = walker.getTableDistances().differences(setA,setB);
    auto & waveValue = walker.getLogWave();

    
    
    const int NA = getN(stateA);
    const int NB=getN(stateB);    
    constexpr int D = getDimensions();
    real_t  tmp,tmp1,tmp2;
    
    int k=0;
    for (int i=0;i<NA;i++)
      {
	for (int j=0;j<NB;j++)
	  {
	    auto d = distances(k);
	    
	    J.evaluateDerivatives(d,tmp,tmp1,tmp2);
	    
	    laplacian+=2*(tmp2 + (D-1)*tmp1/d);
	    waveValue+=tmp;
	    
	    for(int id=0;id<D;id++)
	      {
		gradientA(i,id)+=differences(k,id)/d * tmp1;
		gradientB(j,id)-=differences(k,id)/d * tmp1;
	      }
				
	    k++;
	  }
		
      
      }
  }

  
  static std::string name()   {return "jastrow2bDis/" + jastrow_t::name();}

  virtual std::string print()  override {
   
    return J.print(0,getGeometry().getLBox(0)/2. , 10000)
      ;}

    
private:
  
  int setA;
  int setB;
  jastrow_t J;
  
};

template<class jastrow_t>
class jastrowTwoBodyWavefunctionUndistinguishable : public wavefunction
{
public:
  
  jastrowTwoBodyWavefunctionUndistinguishable(jastrow_t J_,const geometry_t  &geo_, int setA_=0,int setB_=0) : setA(setA_), setB(setB_),J(J_),wavefunction::wavefunction(geo_)
  {
    if (setA != setB) throw invalidInput("setA != setB in undistinguishable two body jastrow");
  }

  jastrowTwoBodyWavefunctionUndistinguishable(const json_t & j,const geometry_t & geo ) : jastrowTwoBodyWavefunctionUndistinguishable( jastrow_t(j["jastrow"]), geo,j["sets"][0].get<int>() ,j["sets"][1].get<int>() ) {}

  virtual std::vector<int> sets() const {return {setA,setB} ;}
  
  
  virtual real_t operator()(const walker_t & walker) 
	{		
	  auto & dis = walker.getTableDistances().distances(setA,setB);
	  
	  int N=size(dis);
	  real_t sum=0;

	  for(int i=0;i<N;i++)
	    {
	      sum+=J.d0(dis(i));
	    }
	  return sum;
	};
  
  virtual void accumulateDerivatives( walker_t & walker ) override
  {
    auto & state = walker.getStates()[setA];

    auto & gradient = walker.getGradients()[setA];
    
    auto & laplacian = walker.getLaplacianLog();
    auto & distances = walker.getTableDistances().distances(setA,setB);
    auto & differences = walker.getTableDistances().differences(setA,setB);
    auto & waveValue = walker.getLogWave();
    
    
    const int N = getN(state);
    constexpr int D = getDimensions();
    
    real_t  tmp,tmp1,tmp2;
    
    int k=0;
    for (int i=0;i<N;i++)
      {
	for (int j=0;j<i;j++)
	  {
	    auto d = distances(k);
	    
	    J.evaluateDerivatives(d,tmp,tmp1,tmp2);
	    
	    laplacian+=2* (tmp2 + (D-1)*tmp1/d);
	    waveValue+=tmp;
	
	    for(int id=0;id<D;id++)
	      {
		gradient(i,id)+=differences(k,id)/d * tmp1;
		gradient(j,id)-=differences(k,id)/d * tmp1;
	      }
				
	    k++;
	  }
		
      
      }
  }

  static std::string name()   {return "jastrow2bUnDis/" + jastrow_t::name();}
  
  
  virtual std::string print()  override {return J.print(0,getGeometry().getLBox(0)/2. , 10000) ;}
  
private:
    
  int setA;
  int setB;
  jastrow_t J;
  
};
