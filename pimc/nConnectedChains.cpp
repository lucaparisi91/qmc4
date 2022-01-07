#include "nConnectedChains.h"
namespace pimc
{
    size_t nConnectedChains::count(const configurations_t & configurations, int iSet )
    {
        const auto & group = configurations.getGroups()[ iSet ];
        
        _visited.resize(group.size() );

        std::fill(_visited.begin(), _visited.end(), false);

        size_t nRings=0;
        
        for ( int i=group.iStart ; i<= group.iEnd ; i++ )
        {
            if ( not _visited[i - group.iStart] )
            {
                int iChain=i;
                nRings+=1;
                do
                {
                    _visited[iChain - group.iStart]=true;
                    iChain=configurations.getChain(iChain).next;
                }
                while ( (iChain>=0) and (iChain!=i) );
            }
        }

        
        return nRings;
    }

};
