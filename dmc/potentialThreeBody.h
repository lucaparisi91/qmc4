#ifndef POTENTIAL_THREE_BODY_H
#define POTENTIAL_THREE_BODY_H


#include "potential.h"
#include "tools3B.h"

template<class functor_t>
class potentialThreeBodyUnDis :  public potential
{
public:
  using potential::operator();
  
  potentialThreeBodyUnDis( const functor_t & V_,int setA_,const geometry_t & geo) : potential(geo),V(V_),setA(setA_){};
  
  potentialThreeBodyUnDis(const json_t & j, const geometry_t & geo) :
    potentialThreeBodyUnDis(functor_t(j["functor"]),j["sets"][0].get<int>(),geo)
  {
    
  }

  
  virtual real_t operator()(const walker_t & w)
  {
    const auto & dis=w.getTableDistances().distances(setA,setA);
    real_t sum=0;

    const int N = getN(w.getStates()[setA]);

    
    LOOP3B( N,
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	    
	   auto R= sqrt(rji*rji + rki*rki + rkj*rkj);
	    
	   sum+=V(R);
	   
	   )

    return sum;
    
  }
  
  static std::string name() {return "potentialThreeBodyUnDis/" + functor_t::name(); }
  
  virtual std::vector<int> sets() const {return {setA,setA,setA};} ;
  
  
private:
  
  functor_t V;
  int setA;
  
};



template<class functor_t>
class potentialThreeBodyDis :  public potential
{
public:
  using potential::operator();
  
  potentialThreeBodyDis( const functor_t & V_,int setA_,int setB_,int setC_,const geometry_t & geo) : potential(geo),V(V_),setA(setA_),setB(setB_),setC(setC_)
  {
     if ( (setA == setB) or (setB == setC) or (setA == setC) )
      {
	throw invalidInput("Sets in three body wavefunction should be distinguishible");
      }
  };

  
  potentialThreeBodyDis(const json_t & j, const geometry_t & geo) :
    potentialThreeBodyDis(functor_t(j["functor"]),j["sets"][0].get<int>(),j["sets"][1].get<int>(),j["sets"][2].get<int>(),geo)
  {
    
  }

  
  virtual real_t operator()(const walker_t & w)
  {
    auto & disAB = w.getTableDistances().distances(setA,setB);
    auto & disBC = w.getTableDistances().distances(setB,setC);
    auto & disAC = w.getTableDistances().distances(setA,setC);   
    
  
    real_t sum=0;
    int NA = getN( w.getStates()[setA] );
    int NB = getN( w.getStates()[setB] );
    int NC = getN( w.getStates()[setC] );
    

    // performs a 3b loop on all particles
    
    LOOP3B_DIS( NA , NB , NC ,

		auto rij = disAB(ij);	
		auto rjk = disBC(jk);	
		auto rik = disAC(ik);
	    
		auto R= std::sqrt(rij*rij + rjk*rjk + rik*rik);
		sum+=V(R);
		) ;
    
    return sum;
    

    
  }
  
  static std::string name() {return "potentialThreeBodyDis/" + functor_t::name(); }
  
  virtual std::vector<int> sets() const {return {setA,setB,setC};} ;
  
  
private:
  
  functor_t V;
  int setA;
  int setB;
  int setC;
  
};








#endif
