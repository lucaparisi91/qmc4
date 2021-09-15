#include "toolsPimc.h"
#include "qmcExceptions.h"

namespace pimc
{

std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads)
    {
        auto [t0,t1] = timeSlice;
        std::array<std::array<int,2>, 2> timeSlices;

        if ( (t1 - t0 ) >= nBeads )
        {
            throw invalidInput("Length of time slice should be samller then the number of beads");
        }

        
        if ( t1 < nBeads )
        {
            if ( t0< 0 )
            {
                timeSlices[0]={nBeads + t0,nBeads-1};
                timeSlices[1]={0,t1};
            }
            else
            {
                timeSlices[0]={t0,t1};
                timeSlices[1]={0,-1};
            }
            
            
        }
        else 
        {
                timeSlices[0]={t0, nBeads-1   };
                timeSlices[1]={0,t1%nBeads};

        }

        return timeSlices;

    }


Real freeParticleLogProbability(std::array<Real,getDimensions()> & delta,Real tau,Real mass)
    {
        const Real D = 0.5;
        Real var = 2 * D * tau / mass;

        Real p=0;
        
        for(int d=0;d<getDimensions();d++)
        {
            p+= delta[d]*delta[d];
        }
        p*=-0.5 /var;
        p+= -0.5*getDimensions()*log(2*M_PI*var);
        
        return p;
    }

Real average(const std::vector<Real> & observables)
{
    Real sum=0;
    for (auto value : observables)
    {
        sum+=value;
    }
    return sum/observables.size();
    
}


int getTotParticlesFromJson(const json_t & j)
{
    auto nParticles=j["particles"].get<std::vector<int> >();

    return std::accumulate(nParticles.begin(),nParticles.end(),0);
}

}