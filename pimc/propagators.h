#include "traits.h"
#include <memory>


namespace pimc
{

    class caoBernePropagator
    {
    public:

        caoBernePropagator(Real timeStep, Real a_) : tau(timeStep),a(a_),D(1) {}
        
        Real logEvaluate( const std::array<Real,3> & x1 , const std::array<Real,3> & x2 ) const
        {
            Real r1 = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
            Real r2 = std::sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
            Real r1Dotr2=x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];


             Real value =  
                        1 - (a*(r1 + r2) - a*a)/(r1*r2) *
                        exp(- 
                            (r1*r2 + a*a -a*(r1+r2) )*(1+r1Dotr2/(r1*r2))/(2*D*tau)
                        );
            

            return -log(value);    

        }


        inline Real logTimeDerivative2( const std::array<Real,3> & x1 , const std::array<Real,3> & x2 ) const
        {
            Real r1 = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
            Real r2 = std::sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
            Real r1Dotr2=x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];

            return -(-a+r1)*(a - r2)*(r1Dotr2 + r1*r2) / (

                    ( -r1*r2 + exp(  (-a+r1)*(-a + r2)*(r1Dotr2 + r1*r2)/(2*tau*D*r1*r2)) * r1*r1*r2*r2/(a*(-a + r1 + r2))
                    ) * 
                    (2*D*tau*tau)
            )   ; 
        }

        inline Real logTimeDerivative( const std::array<Real,3> & x1 , const std::array<Real,3> & x2 ) const
        {
            Real r1 = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
            Real r2 = std::sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
            Real r1Dotr2=x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
            Real cosTeta = r1Dotr2/(r1*r2);

            Real E=exp(-(r1*r2 + a*a - a*(r1+r2) )*(1+cosTeta)/(2*D*tau));


            Real dgdtau= -a*(r1 + r2 - a)/(r1*r2) *(r1*r2 + a*a -a*(r1+r2)) *(1+cosTeta)/(2*D*tau*tau) * E;

            
            Real g=1 - a*( r1 + r2 - a)/(r1*r2)*E;

            return -dgdtau/g;

        }

        Real logGradientLeft2(const std::array<Real,3> & x1, const std::array<Real,3> & x2, int d) const
        {
            Real r1 = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
            Real r2 = std::sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
            Real r1Dotr2=x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];

            return 
            (
            -exp( (-a+r1)*(a- r2) *(r1Dotr2 + r1*r2)/(2*D*tau*r1*r2)    )*(

                -a*x1[d]/(r1*r1) + a * x1[d]*(-a + r1 + r2)/(r1*r1*r1)  +

                a*(-a + r1 + r2)*(-a + r2)*
                (  -a*x2[d]*(r1*r1) + x1[d]*a*r1Dotr2 + r1*r1*(x2[d]*r1 + x1[d]*r2)  
            ) / (2 * D * tau * r1* r1*r1*r1*r2) ) 
            )/
            (
                r2*(1 - 
                    a *exp( (-a + r1)*(a-r2)*(r1Dotr2 + r1*r2)/(2*D*tau*r1*r2)  )*(-a + r1 + r2)/(r1*r2)
                )
            );
        }

        Real logGradientLeft(const std::array<Real,3> & x1, const std::array<Real,3> & x2, int d) const
        {
            Real r1 = std::sqrt(x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2]);
            Real r2 = std::sqrt(x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2]);
            Real r1Dotr2=x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];

            Real cosTeta = r1Dotr2/(r1*r2);

            Real E=exp(-(r1*r2 + a*a - a*(r1+r2) )*(1+cosTeta)/(2*D*tau));

            Real dgdcosTeta = a*(r1 + r2 - a)/(r1*r2) * (r1*r2 + a*a - a *( r1 + r2)  )/(2*D*tau) * E;
            Real dgdr = -a / (r2*r1) *( 1 - (r1 + r2 - a)/r1 + (r1 + r2 -a)*(a - r2)*(1+ cosTeta)/(2*D*tau))*E;


            Real drdrd = x1[d]/r1;
            Real dcosTetadrd = ( x2[d]/(r1*r2) - x1[d]*r1Dotr2/(r2*r1*r1*r1) );

            Real dgdrd = dgdr * drdrd + dgdcosTeta*dcosTetadrd;

            Real g=1 - a*( r1 + r2 - a)/(r1*r2)*E;
            
            return -dgdrd/g;
            
        }



         Real logGradientRight(const std::array<Real,3> & x1, const std::array<Real,3> & x2, int d) const
        {
            return logGradientLeft(x2,x1,d);
        }

    private:

    Real D;
    Real tau;
    Real a;
    };
}