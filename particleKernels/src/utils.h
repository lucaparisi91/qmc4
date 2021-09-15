#ifndef UTILS_H
#define UTILS_H

using Real = double;
#include <array>
#include <random>

namespace particleKernels
{
namespace utils
{
     inline Real differencePBC(Real t, Real lBox, Real lBoxInverse ) {return ( t - std::round(t*lBoxInverse )*lBox);}

     
    Real restrictToBox(Real x, Real left, Real lBox, Real lBoxInverse);

    void restrictToBox(Real * positionsPBC, const Real * positionsOld, int iStart,int iEnd, int N, const std::array<Real,3> & left, const std::array<Real,3> & lBox );
    
    void initRandom(Real * particles,int dims,int iStart,int iEnd,size_t N,  std::default_random_engine & generator,
 Real xmin=0,Real xmax=1);

 
 template<class V_t,int dims >
Real twoBodyDistancesIsotropicReduction(const __restrict Real * particles,const V_t & op,size_t N,const std::array<Real,dims> &  lBox)
{
    Real sum2b=0;

    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(runtime)
    for (int iParticle=0;iParticle<N;iParticle++)
    {   
        for (int jParticle=0;jParticle<iParticle;jParticle++)
        {
            Real r2=0;
            for(int d=0;d<dims;d++)
            {
                Real diffd=( particles[iParticle + d*N ] - particles[jParticle + d*N] );
                diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                r2+= diffd * diffd;
            }


            sum2b+=op(std::sqrt(r2));

        }
    }

return sum2b;

}



struct gaussianInteractionWithCutOff
{
    gaussianInteractionWithCutOff(Real alpha, Real cutOff) :
    _alpha(alpha),
    _cutOff(cutOff)
    {}
    inline Real operator()(Real r) const
    {
        if (r < _cutOff)
        {
            return std::exp(-_alpha*r*r);
        }
        else
        {
            return 0;
        }
        
    }

    private:
    Real _alpha;
    Real _cutOff;
};


struct tetaInteraction
{
    tetaInteraction(Real V0,Real cutOff) :
    _cutOff(cutOff),
    _V0(V0) {}

    inline Real operator()(Real r) const
    {
        if (r < _cutOff)
        {
            return _V0;
        }
        else
        {
            return 0;
        }
        
    }

    private:
    Real _V0;
    Real _cutOff;
};

};
}

#endif