#include "forces.h"


namespace pimc
{
    void totalForce::init(int M, int nChains)
    {
        forceData = std::make_shared<tensor_t>(nChains,getDimensions(),M+1);
        
    }

    void totalForce::setConstant(Real alpha,const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange)
    {
        auto & force = *forceData;


        
         for (int t=timeRange[0] ;t<= timeRange[1] ; t++  )
                for (int i=particleRange[0];i<=particleRange[1];i++ )
                    {
                        for (int d=0;d<getDimensions();d++)
                        {
                            force(i,d,t)=0;
                        }
                    }                   
    }


    void totalForce::evaluate(const pimcConfigurations_t & configurations,const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange)
    {
        setConstant(0,timeRange,particleRange);

        for (const auto & currentForce : forces)
        {
            currentForce->add(configurations,timeRange,particleRange,*forceData);
        }
    }


}