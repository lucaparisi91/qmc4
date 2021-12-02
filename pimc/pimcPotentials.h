#ifndef PIMCPOTENTIALS_H
#define PIMCPOTENTIALS_H

#include "toolsPimc.h"
#include <cmath>
#include "qmcExceptions.h"
namespace pimc
{

using Real = double;

class poschlTellerPotential
{
    public:

    poschlTellerPotential(const json_t & j) : 
    poschlTellerPotential::poschlTellerPotential(j["V0"].get<Real>(),j["R"].get<Real>()  ) {}


    poschlTellerPotential(Real V0_,Real R_) : V0(V0_),R(R_) {}
    
    #if DIMENSIONS == 3
    Real operator()(Real x,Real y, Real z)  const { Real r = std::sqrt(x*x + y*y + z*z); return  -V0/std::pow(cosh(r/R),2);}

    Real gradX(Real x,Real y , Real z) const {Real r = std::sqrt(x*x + y*y + z*z); return x/r *radialDerivative(r) ;   }

    Real gradY(Real x,Real y , Real z) const {Real r = std::sqrt(x*x + y*y + z*z); return y/r *radialDerivative(r) ;   }

    Real gradZ(Real x,Real y , Real z) const {Real r = std::sqrt(x*x + y*y + z*z); return z/r *radialDerivative(r) ;   }

    #endif

    #if DIMENSIONS == 2
    Real operator()(Real x,Real y) const { Real r = std::sqrt(x*x + y*y); return  -V0/std::pow(cosh(r/R),2);}

    Real gradX(Real x,Real y) const {Real r = std::sqrt(x*x + y*y ); return x/r *radialDerivative(r) ;   }

    Real gradY(Real x,Real y ) const {Real r = std::sqrt(x*x + y*y ); return y/r *radialDerivative(r) ;   }

    #endif

    #if DIMENSIONS == 1

    Real operator()(Real x) const { Real r = std::abs(x); return  -V0/std::pow(cosh(r/R),2);}

    Real gradX(Real x) const {Real r = std::abs(x ); return x/r *radialDerivative(r) ;   }


    #endif




    private:

    Real radialDerivative(Real r) const { Real x=r/R; return 2*V0/R*sinh(x)/(std::pow(cosh(x),3)) ;}


    Real V0;
    Real R;

};


class harmonicPotential
{
    public:
    

    harmonicPotential(const json_t & j) :
    harmonicPotential::harmonicPotential(j["omega"].get<std::vector<Real> >() ) {}

    static  std::string name() {return "harmonic";}


    harmonicPotential(std::vector<Real> omegas_={0.5,0.5,0.5}) : 
 omega(omegas_) {
     if ( omega.size() != getDimensions() )
     {
         throw invalidInput("Number of dimenions of omega vector is inconsistent");
     }
 }


    #if DIMENSIONS == 3

    Real operator()(Real x,Real y,Real z) const { return omega[0]*x*x + omega[1]*y*y + omega[2]*z*z;}

    Real gradX(Real x,Real y,Real z) const {return 2*omega[0] *x   ;}
    Real gradY(Real x,Real y,Real z) const {return 2*omega[1] *y   ;}
    Real gradZ(Real x,Real y,Real z) const {return 2*omega[2] *z   ;}

    #endif


    #if DIMENSIONS == 2
    Real operator()(Real x,Real y) { return omega[0]*x*x + omega[1]*y*y;}

    Real gradX(Real x,Real y) {return 2*omega[0] *x   ;}
    Real gradY(Real x,Real y) {return 2*omega[1] *y   ;}

    #endif


    #if DIMENSIONS == 1
    Real operator()(Real x) const { return omega[0]*x*x;}
    Real gradX(Real x) const {return 2*omega[0] *x   ;} 

    #endif

    private:


    std::vector<Real> omega;

};

class isotropicHarmonicPotential
{
    public:
    
    isotropicHarmonicPotential(const json_t & j) :
    isotropicHarmonicPotential(j["omega"].get<Real>() ) {}


    static  std::string name() {return "isotropicHarmonic";}


    isotropicHarmonicPotential(Real omega) : 
 _omega(omega) {

 }

    Real radialDerivative(Real r) const {return _omega*r;}

    Real operator()(Real r) const { return 0.5*_omega*r*r;}
    

    #if DIMENSIONS == 3


    Real operator()(Real x,Real y,Real z) const { return 0.5*(x*x + y*y + z*z) ;}

    Real gradX(Real x,Real y,Real z) const {return _omega*x   ;}
    Real gradY(Real x,Real y,Real z) const {return _omega*y   ;}
    Real gradZ(Real x,Real y,Real z) const {return _omega*z   ;}
    #endif

    #if DIMENSIONS == 1
    Real gradX(Real x) const {return _omega*x   ;}
    #endif



    private:

    Real _omega;
};

class gaussianPotential
{
    public:

    gaussianPotential(const json_t & j) :
    gaussianPotential::gaussianPotential(j["V0"].get<Real>() , j["alpha"].get<Real>() ) {}

    static  std::string name() {return "gaussian";}


    gaussianPotential(  Real V0, Real alpha) 
    :   _V0(V0),_alpha(alpha) {}

    Real radialDerivative(Real r) const {return -2*_alpha*r*_V0*exp(-_alpha*r*r);}
   


    #if DIMENSIONS == 2
    Real operator()(Real x, Real y) const { Real r2= x*x + y*y; return _V0*exp(-_alpha*r2);}
    Real gradX(Real x,Real y) { Real r2= x*x + y*y; return -2*x*_V0*_alpha*exp(-_alpha*r2);}
    Real gradY(Real x,Real y) { Real r2= x*x + y*y; return -2*y*_V0*_alpha*exp(-_alpha*r2);}
    #endif


    #if DIMENSIONS == 1
    Real operator()(Real x) const { return _V0*exp(-_alpha*x*x);}
    Real gradX(Real x) const {return -2*x*_V0*_alpha*exp(-_alpha*x*x);}
    #endif



    #if DIMENSIONS == 3

    Real operator()(Real r) const {
         Real r2= r*r; 
         return _V0*exp(-_alpha*r2);}


    Real operator()(Real x, Real y,Real z) const { Real r2= x*x + y*y + z*z; return _V0*exp(-_alpha*r2);}
    Real gradX(Real x,Real y,Real z) const { Real r2= x*x + y*y + z*z; return -2*x*_V0*_alpha*exp(-_alpha*r2);}
    Real gradY(Real x,Real y,Real z) const { Real r2= x*x + y*y + z*z; return -2*y*_V0*_alpha*exp(-_alpha*r2);}
    Real gradZ(Real x,Real y,Real z) const { Real r2= x*x + y*y + z*z; return -2*z*_V0*_alpha*exp(-_alpha*r2);}
    #endif

    private:

    Real _alpha;
    Real _V0;

};

}

#endif