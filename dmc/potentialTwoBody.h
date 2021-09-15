#ifndef POTENTIAL_TWO_BODY_H
#define POTENTIAL_TWO_BODY_H


#include "potential.h"

template<class functor_t>
class potentialTwoBody :  public potential
{
public:
  using potential::operator();
  
  potentialTwoBody( const functor_t & V_,int setA_,int setB_,const geometry_t & geo) : potential(geo),V(V_),setA(setA_),setB(setB_){};

  
  potentialTwoBody(const json_t & j, const geometry_t & geo) :
    potentialTwoBody(functor_t(j["functor"]),j["sets"][0].get<int>(),j["sets"][1].get<int>(),geo)
  {
    
  }

  
  virtual real_t operator()(const walker_t & w)
  {
    const auto & dis=w.getTableDistances().distances(setA,setB);
    real_t sum=0;
    for (int i=0;i<dis.size();i++)
      {
	sum+=V(dis(i));
      }
    return sum;
    
  }
  
  static std::string name() {return "potentialTwoBody/" + functor_t::name(); }
  
  
  virtual std::vector<int> sets() const {return {setA,setB};} ;
  
  
private:
  
  functor_t V;
  int setA;
  int setB;
  
};

#endif
