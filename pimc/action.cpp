#include "action.h"

namespace pimc
{

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeRange , std::array<int,2> particleRange  )
{
    const auto & data = configurations.dataTensor();
    auto & geo=getGeometry();

    Real sum=0;

    assert(timeRange[0]>=0);
    assert(timeRange[1]<configurations.nBeads() );

    assert(particleRange[0] >=0  );
    assert(particleRange[1] <= configurations.nChains() );
#ifndef DEBUG
#endif

   
    for (int i=particleRange[0];i<=particleRange[1];i++ )
        for(int t=timeRange[0];t<=timeRange[1] ; t++ )
            {
                for(int d=0;d<getDimensions();d++)
                {
                    sum+=
                    (data(i,d,t+1 ) - data(i,d,t))*(data(i,d,t+1 ) - data(i,d,t));
                    
                    //geo.difference(data(i,d,t+1 ) - data(i,d,t),d) *
                    //geo.difference(data(i,d,t+1 ) - data(i,d,t),d);
                }
                
            }
            
    return sum/(4* D * getTimeStep());
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChange )
{   
    Real sum=0;

    const auto & particleGroups = configurations.getGroups();

    for (const auto & group : particleGroups)
    {
        if ( ( iChange >= group.iStart )and iChange <= group.iEnd  )
        {
            sum+=group.mass*evaluate(configurations,timeSlices, {iChange,iChange});
        }
    }
    return sum;
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations , std::array<int,2> timeSlices , int iChain1 , int iChain2 )
{   
    return evaluate(configurations,timeSlices,iChain1) + evaluate(configurations,timeSlices,iChain2) ; 
}

Real kineticAction::evaluate( pimcConfigurations_t & configurations )
{
    Real sum=0;
    
    const auto & particleGroups = configurations.getGroups();

    for (auto & group : particleGroups)
    {
        sum+=evaluate(configurations  , {0 , nBeads-1} , {group.iStart,group.iEnd} ) ;
    }

    return sum;
};


}