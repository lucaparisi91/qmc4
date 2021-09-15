#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "traits.h"
#include <tuple>
#include <cmath>
#include <initializer_list>


class geometry
{
public:

	geometry(){}
	using particles_t = ::state_t;
	using diff_t =  ::difference_t;
	
	virtual diff_t differencesTwoBody(const particles_t & data) const = 0; // distances between undistinguishible particles

	virtual diff_t differencesTwoBody(const particles_t & data1, const particles_t & data2) const = 0; // distances between distinguishiblle particles

	virtual diff_t differencesOneBody(const particles_t & data , const std::array<real_t,3> & x) const = 0;

	virtual diff_t differencesOneBody(const particles_t & data, std::initializer_list<real_t> l) const;

	virtual diff_t differences(const particles_t & data, std::initializer_list<real_t> l) const { return differencesOneBody(data,l);}

	virtual diff_t differences(const particles_t & data) const { return differencesTwoBody(data);}

	virtual diff_t differences(const particles_t & data1, const particles_t & data2 ) const { return differencesTwoBody(data1,data2);}
        virtual real_t getLBox(int i) const {return 0;}
};


real_t norm(real_t x,real_t y, real_t z);

distance_t norm( const difference_t & diffs);


class geometryPBC : public geometry
{
public:

  geometryPBC(){}
  geometryPBC(real_t lBoxx_,real_t lBoxy_,real_t lBoxz_) : lBox{ lBoxx_,lBoxy_,lBoxz_},lBoxInverse{1./lBoxx_,1./lBoxy_,1./lBoxz_} {} ;


  geometryPBC(real_t lBox ) : geometryPBC(lBox,lBox,lBox) {};
  
  real_t difference(real_t t, int i ) const {return ( t - std::round(t*lBoxInverse[i] )*lBox[i]);}

  auto difference(real_t diffx,real_t diffy,real_t diffz) const { return std::make_tuple( difference(diffx,0) , difference(diffy,1) , difference(diffz,1)   );}

  virtual diff_t differencesOneBody(const particles_t & particleData, const std::array<real_t,3> & x) const;

  virtual diff_t differencesTwoBody(const particles_t & particleData) const ;

  virtual diff_t differencesTwoBody(const particles_t & data1,const particles_t & data2) const ;

  virtual real_t getLBox(int i) const override {return lBox[i];}

  virtual void differencesTwoBody(diff_t & differences , const particles_t & particleData, int i) const; // updates differences when particle i is changed

private:
	real_t lBox [3];
	real_t lBoxInverse [3];

};

#endif
