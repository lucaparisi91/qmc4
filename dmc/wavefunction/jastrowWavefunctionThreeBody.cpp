#include "walkers.h"
#include "jastrowWavefunctionThreeBody.h"



hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(int setA_,real_t R0) : setA(setA_),V(R0) {}


hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(const json_t & j) :
  hardSphereConstraintThreeBodyUnDis::hardSphereConstraintThreeBodyUnDis(j["sets"][0].get<real_t>(), j["hardRadius"].get<real_t>() )
{
  
}

bool hardSphereConstraintThreeBodyUnDis::operator()(const walker_t & w)
  {

    const auto & dis=w.getTableDistances().distances(setA,setA);
    
    const int N = getN(w.getStates()[setA]);
    
    LOOP3B( N,
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	    
	   auto R= sqrt(rji*rji + rki*rki + rkj*rkj);
	    
	    if (R <= V.R0() )
	      {
		return false;
	      }
	    )
      
    return true;  
  }



hardSphereConstraintThreeBodyDis::hardSphereConstraintThreeBodyDis(int setA_,int setB_,int setC_,real_t R0) : setA(setA_),setB(setB_),setC(setC_),V(R0) {}


hardSphereConstraintThreeBodyDis::hardSphereConstraintThreeBodyDis(const json_t & j) :
  hardSphereConstraintThreeBodyDis::hardSphereConstraintThreeBodyDis(j["sets"][0].get<real_t>(),j["sets"][1].get<real_t>(),j["sets"][2].get<real_t>(), j["hardRadius"].get<real_t>() )
{
  
}


bool hardSphereConstraintThreeBodyDis::operator()(const walker_t & w)
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

		auto R= sqrt(rij*rij + rjk*rjk + rik*rik);
		
		if (R <= V.R0() )
		  {
		    return false;
		  }
		)
      
      return true;  
  }

