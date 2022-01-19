
#ifndef PRIMITIVEPROPAGATOR_H
#define PRIMITIVEPROPAGATOR_H

#include "traits.h"

namespace pimc
{
    template<class V_t>
    class primitivePropagator
    {
    public:

        primitivePropagator(Real timeStep, V_t V) : tau(timeStep)
        {
            _V=std::make_unique<V_t>(V);

        }

        primitivePropagator(const json_t & j) {
            throw std::runtime_error("Constructor from json file not implemented");
        }

        

        Real logEvaluate( const std::array<Real,DIMENSIONS> & x1 , const std::array<Real,DIMENSIONS > & x2 ) const
        {
            Real r2_1=0,r2_2=0;
            for(int d=0;d<DIMENSIONS;d++)
            {
                r2_1+=x1[d]*x1[d];
                r2_2+=x2[d]*x2[d];
            }
            Real r1=std::sqrt(r2_1);
            Real r2=std::sqrt(r2_2);
            
            return tau*0.5*( (*_V)(r1) + (*_V)(r2) );
        }

        inline Real logTimeDerivative( const std::array<Real,DIMENSIONS> & x1 , const std::array<Real,DIMENSIONS> & x2 ) const
        {
            Real r2_1=0,r2_2=0;
            for(int d=0;d<DIMENSIONS;d++)
            {
                r2_1+=x1[d]*x1[d];
                r2_2+=x2[d]*x2[d];
            }
            Real r1=std::sqrt(r2_1);
            Real r2=std::sqrt(r2_2);

            return 0.5 * ( (*_V)(r1) + (*_V)(r2) );
        }
        
        Real logGradientLeft(const std::array<Real,DIMENSIONS> & x1, const std::array<Real,DIMENSIONS> & x2, int d) const
        {
            Real r2_1=0;
            for(int d=0;d<DIMENSIONS;d++)
            {
                r2_1+=x1[d]*x1[d];
                
            }
            Real r1=std::sqrt(r2_1);
            
            return 0.5 * tau *  (_V->radialDerivative(r1) ) * x1[d]/r1; 
        }


        Real logGradientRight(const std::array<Real,DIMENSIONS> & x1, const std::array<Real,DIMENSIONS> & x2, int d) const
        {
            return logGradientLeft(x2,x1,d);
        }


    private:

    std::unique_ptr<V_t> _V;
    Real tau;
    };
}

#endif