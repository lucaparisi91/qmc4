#include "actionTwoBody.h"
namespace pimc
{


Real actionTwoBody::evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange, int iParticle)
    {
        auto const & groupA = configurations.getGroups()[setA];
        auto const & groupB = configurations.getGroups()[setB];

        // if modified particle is not involved in the evaluate of the action
        bool isInA = groupA.contains(iParticle);
        bool isInB = groupB.contains(iParticle);

        if (not (isInA or isInB) ) {return 0;} 
        if (setA == setB)
        {
            return _kernel->evaluateTriangular(configurations.dataTensor() , timeRange,{iParticle,iParticle} , {groupA.iStart,groupA.iEnd} );
        }
        else
        {
            std::array<int,2> rangeB;

            if (isInA)
            {
                rangeB = std::array<int,2>{groupB.iStart,groupB.iEnd };
            }
            else
            {
                rangeB = std::array<int,2>{groupA.iStart,groupA.iEnd };
            }

            return _kernel->evaluateRectangular(configurations.dataTensor(), timeRange , {iParticle,iParticle} , rangeB ) ;
        }
    }


Real actionTwoBody::evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange)
    {
        auto const & groupA = configurations.getGroups()[setA];
        auto const & groupB = configurations.getGroups()[setB];

        auto rangeA = std::array<int,2>{groupB.iStart,groupB.iEnd };

        if (setA == setB)
        {
            return _kernel->evaluateTriangular(configurations.dataTensor() , timeRange,{groupA.iStart,groupA.iEnd} , {groupA.iStart,groupA.iEnd} );
        }
        else
        {
            auto rangeB = std::array<int,2>{groupB.iStart,groupB.iEnd };

            return _kernel->evaluateRectangular(configurations.dataTensor(), timeRange , rangeA , rangeB ) ;
        }
    }

}