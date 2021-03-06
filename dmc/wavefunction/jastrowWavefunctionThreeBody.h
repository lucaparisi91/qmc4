#ifndef JASTROWWAVEFUNCTIONTHREEBODY_H
#define JASTROWWAVEFUNCTIONTHREEBODY_H

#include "wavefunction/wavefunction.h"

#include "tools3B.h"
#include "hardSphere.h"

class hardSphereConstraintThreeBodyUnDis : public constraint

{
public:
  hardSphereConstraintThreeBodyUnDis(int setA_,real_t R0);

  hardSphereConstraintThreeBodyUnDis(const json_t & j) ;
  
  virtual bool operator()(const walker_t & w) override;
private:
  int setA;
  hardSphere V;
};


class hardSphereConstraintThreeBodyDis : public constraint

{
public:
  hardSphereConstraintThreeBodyDis(int setA_,int setB_,int setC_,real_t R0);

  hardSphereConstraintThreeBodyDis(const json_t & j) ;
  
  virtual bool operator()(const walker_t & w) override;
  
private:
  int setA,setB,setC;
  hardSphere V;
};



template<class jastrow_t>
class jastrowThreeBodyWavefunctionUnDistinguishable : public wavefunction
{
public:
  
  jastrowThreeBodyWavefunctionUnDistinguishable(jastrow_t J_,const geometry_t  &geo_, int setA_=0) : setA(setA_),J(J_),wavefunction::wavefunction(geo_)
  {
    
  }

  
  jastrowThreeBodyWavefunctionUnDistinguishable(const json_t & j,const geometry_t & geo ) : jastrowThreeBodyWavefunctionUnDistinguishable( jastrow_t(j["jastrow"]), geo,j["sets"][0].get<int>()  )
  {
    if ( j.contains("hardSphereRadius") )
	   {
	     auto R0 = j["hardSphereRadius"].get<real_t>();
	     
	     addHardSphereConstraint( R0 );
	   }
      
  }

  
  virtual real_t operator()(const walker_t & walker)
  {
    auto & dis = walker.getTableDistances().distances(setA,setA);  
  
    real_t sum=0;
    int N = getN( walker.getStates()[setA] );

    // performs a 3b loop on all particles
    
    LOOP3B( N ,
	    auto rji = dis(ji);	
	    auto rkj = dis(kj);	
	    auto rki = dis(ki);
	    
	    auto R= std::sqrt(rji*rji + rki*rki + rkj*rkj);
	    sum+=J.d0(R);
	    ) ;
  
    return sum;
    
  }


  
  
  virtual std::vector<int> sets() const {return {setA,setA} ;}
  
  virtual void accumulateDerivatives( walker_t & walker ) override
  {
    auto & state = walker.getStates()[setA];
    
    auto & gradient = walker.getGradients()[setA];
    
    auto & laplacian = walker.getLaplacianLog();
    auto & dis = walker.getTableDistances().distances(setA,setA);
    auto & diff = walker.getTableDistances().differences(setA,setA);
    auto & waveValue = walker.getLogWave();
    
    
    const int N = getN(state);
    constexpr int D = getDimensions();
    
    real_t  d0,d1,d2;
    
    
    LOOP3B(N,

	   
	   auto rji = dis(ji);	
	   auto rkj = dis(kj);	
	   auto rki = dis(ki);
	   
	   auto R= std::sqrt(rji*rji + rki*rki + rkj*rkj);
	   
	   J.evaluateDerivatives(R,d0,d1,d2);
	   
	   
	   for (int id=0;id < D;id++)
	     {    
	       gradient(i,id)-=  d1*( diff(ji,id) + diff(ki,id))/R;
	       gradient(j,id)+=  d1*(   diff(ji,id) - diff(kj,id) )/R;
	       gradient(k,id)+=  d1*(   diff(ki,id) + diff(kj,id) )/R;
	     }
	   
	   laplacian+= d2 * 3 ;
	   laplacian += 3 * d1 / R *  (2*D-1);
	   
	   waveValue+=d0;
	   
	   
	   

	   )
  }
  
  static std::string name()   {return "jastrow3bUnDis/" + jastrow_t::name();}
  
  virtual std::string print()  override {
    
    return J.print(0,getGeometry().getLBox(0)/2. , 100000)
      ;}  

  
  void addHardSphereConstraint(real_t R0)
  {
    this->addConstraint(  std::make_unique<hardSphereConstraintThreeBodyUnDis>(  hardSphereConstraintThreeBodyUnDis(setA,R0) ) );
  }
  
private:
  int setA ;
  jastrow_t J;
  
};


template<class jastrow_t>
class jastrowThreeBodyWavefunctionDistinguishable : public wavefunction
{
public:
  
  jastrowThreeBodyWavefunctionDistinguishable(jastrow_t J_,const geometry_t  &geo_, int setA_,int setB_,int setC_) : setA(setA_),setB(setB_),setC(setC_),J(J_),wavefunction::wavefunction(geo_)
  {
    
    if ( (setA == setB) or (setB == setC) or (setA == setC) )
      {
	throw invalidInput("Sets in three body wavefunction should be distinguishible");
      }


  }

  
  jastrowThreeBodyWavefunctionDistinguishable(const json_t & j,const geometry_t & geo ) : jastrowThreeBodyWavefunctionDistinguishable( jastrow_t(j["jastrow"]), geo,j["sets"][0].get<int>() ,j["sets"][1].get<int>() , j["sets"][2].get<int>() )
  {
        if ( j.contains("hardSphereRadius") )
	   {
	     auto R0 = j["hardSphereRadius"].get<real_t>();
	     
	     addHardSphereConstraint( R0 );
	     
	   }
	
  }

  
  virtual real_t operator()(const walker_t & walker)
  {
    auto & disAB = walker.getTableDistances().distances(setA,setB);
    auto & disBC = walker.getTableDistances().distances(setB,setC);
    auto & disAC = walker.getTableDistances().distances(setA,setC);
    
    
  
    real_t sum=0;
    int NA = getN( walker.getStates()[setA] );
    int NB = getN( walker.getStates()[setB] );
    int NC = getN( walker.getStates()[setC] );
    

    // performs a 3b loop on all particles
    
    LOOP3B_DIS( NA , NB , NC ,

		auto rij = disAB(ij);	
		auto rjk = disBC(jk);	
		auto rik = disAC(ik);
	    
		auto R= std::sqrt(rij*rij + rjk*rjk + rik*rik);
		sum+=J.d0(R);
		) ;
    
    return sum;
    
  }

  
  virtual std::vector<int> sets() const {return {setA,setB,setC} ;}
  
  virtual void accumulateDerivatives( walker_t & walker ) override
  {
    auto & stateA = walker.getStates()[setA];
    auto & stateB = walker.getStates()[setB];
    auto & stateC = walker.getStates()[setC];
    
    auto & gradientA = walker.getGradients()[setA];
    auto & gradientB = walker.getGradients()[setB];
    auto & gradientC = walker.getGradients()[setC];

    
    
    
    auto & laplacian = walker.getLaplacianLog();

    auto & disAB = walker.getTableDistances().distances(setA,setB);
    auto & disBC = walker.getTableDistances().distances(setB,setC);
    auto & disAC = walker.getTableDistances().distances(setA,setC);
    
    
    auto & diffAB = walker.getTableDistances().differences(setA,setB);
    auto & diffBC = walker.getTableDistances().differences(setB,setC);
    auto & diffAC = walker.getTableDistances().differences(setA,setC);
    
    auto & waveValue = walker.getLogWave();
    
    
    const int NA = getN(stateA);
    const int NB = getN(stateB);
    const int NC = getN(stateC);
    

    
    constexpr int D = getDimensions();
    
    real_t  d0,d1,d2;
    
    LOOP3B_DIS(NA, NB , NC ,
	       
	       auto rij = disAB(ij);	
	       auto rjk = disBC(jk);	
	       auto rik = disAC(ik);
	       
	       auto R= std::sqrt(rij*rij + rjk*rjk + rik*rik);
	       
	       J.evaluateDerivatives(R,d0,d1,d2);
	   
	       for (int id=0;id < D;id++)
		 {    
		   gradientA(i,id)+=  d1*( diffAB(ij,id) + diffAC(ik,id))/R;
		   gradientB(j,id)+=  d1*(  - diffAB(ij,id) + diffBC(jk,id) )/R;
		   gradientC(k,id)-=  d1*(   diffAC(ik,id) + diffBC(jk,id) )/R;
		 }
	   
	       laplacian+= d2 * 3 ;
	       laplacian += 3 * d1 / R *  (2*D-1);
	       
	       waveValue+=d0;

	   )
  }
  
  static std::string name()   {return "jastrow3bDis/" + jastrow_t::name();}
  
  virtual std::string print()  override {
    
    return J.print(0,getGeometry().getLBox(0)/2. , 100000)
      ;}

  
  void addHardSphereConstraint(real_t R0)
  {
    this->addConstraint(  std::make_unique<hardSphereConstraintThreeBodyDis>(  hardSphereConstraintThreeBodyDis(setA,setB,setC,R0) ) );
  }


  
private:
  int setA ;
  int setB;
  int setC;
  jastrow_t J;
  
};



#endif
