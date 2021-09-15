#ifndef POTENTIAL_H
#define POTENTIAL_H

#include "traits.h"
#include "qmcExceptions.h"
#include "walkers.h"

class walker;

class tableDistances;

class potential
{
	/*
	Represents the logarithm of the wavefunction. Supports wavefunctions acting on at most two different sets
	*/
public:
  using walker_t = walker;
  using states_t=::states_t;
  using state_t = ::state_t;

  potential(const geometry_t & geo_ );

  virtual real_t operator()(const walker_t & state ) = 0;
  const auto & getGeometry() {return *geo;}
  virtual std::vector<int> sets() const = 0 ;

  virtual bool satisfyConstrains(const walker_t & w) {return true;}

  
private:
  const geometry_t * geo;
};




class sumPotentials
{
public:
  using walker_t = potential::walker_t;
  sumPotentials(std::vector<potential*> potentials);
  const auto & potentials() const {return _potentials;}
  virtual real_t operator()(const walker_t & state );
  size_t size(){return _potentials.size();}
private:
  std::vector<potential*> _potentials;
};


class potential1b : public potential
{
public:

  potential1b(const geometry_t & geo, int setA_);
  
  virtual real_t operator()(const walker_t & state)=0;
  
  virtual std::vector<int> sets() const {return {_setA};}

  int setA() const {return _setA;}
  
private:
	int _setA;
};

class harmonicPotential :  public potential1b
{
public:
  using potential1b::operator();
  
  harmonicPotential(const geometry_t & geo,real_t freq , int setA ) ;
  
  harmonicPotential(const json_t & j, const geometry_t & geo);
  
  virtual real_t operator()(const walker_t & state) ;
  
  static std::string name() {return "harmonic1b";}
private:

	real_t omega;
};

class emptyPotential :  public potential1b
{
public:
  using potential1b::operator();

  emptyPotential(const geometry_t & geo) : potential1b(geo,0) {} ;

  virtual real_t operator()(const walker_t & w)  {return 0.;}
};

class squareWellPotential2b :  public potential
{
public:
  using potential::operator();
  
  squareWellPotential2b(const geometry_t & geo,real_t V0_ , real_t R0_ , int setA_, int setB_ ) ;
  
  squareWellPotential2b(const json_t & j, const geometry_t & geo);
  
  virtual real_t operator()(const walker_t & state) ;
  
  static std::string name() {return "squareWell2b";}

  virtual std::vector<int> sets() const {return {setA,setB};} ;
private:
  real_t R0;
  real_t V0;
  int setA;
  int setB;
};


class dipolarPotential2b :  public potential
{
public:
  using potential::operator();
  
  dipolarPotential2b(const geometry_t & geo,real_t D, int setA_, int setB_ ) ;
  
  dipolarPotential2b(const json_t & j, const geometry_t & geo);
  
  virtual real_t operator()(const walker_t & state) ;
  
  static std::string name() {return "dipolar2b";}

  virtual std::vector<int> sets() const {return {setA,setB};} ;
private:
  real_t D;
  int setA;
  int setB;
  int k;  
  
};

class poschTellerPotential2b :  public potential
{
public:
  using potential::operator();
  
  poschTellerPotential2b(const geometry_t & geo,real_t R0, int setA_, int setB_ ) ;
  
  poschTellerPotential2b(const json_t & j, const geometry_t & geo);
  
  virtual real_t operator()(const walker_t & state) ;
  
  static std::string name() {return "poschTeller2b";}
  

  virtual std::vector<int> sets() const {return {setA,setB};} ;
private:
  real_t R0;
  int setA;
  int setB;
  real_t alpha;
  
};








#endif
