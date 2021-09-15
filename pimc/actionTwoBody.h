#include "memory"


class actionTwoBody
{
    public:


    evaluate(const configurations_t & configurations,std::array<int,2> timeRange, const std::array<int,2> & timeRange , int iParticle)
    {
        const groupA = configurations.groups()[setA];
        const groupB = configurations.groups()[setB];


        // if modified particle is not involved in the evaluate of the action
        bool isInA = groupA.contains(iParticle);
        bool isInB = groupB.contains(iParticle);
        if (not (isInA or isInB) ) {return 0;} 
        if (setA == setB)
        {
            
            _kernel->evaluateTriangular(configurations.data() , timeRange,{iParticle,iParticle} , {groupA.iStart,groupA.iEnd} );
        }
        else
        {
            std::array<int,2> range = isInB ? {groupA.iStart,groupA.iEnd} :{groupB.iStart,groupB.iEnd} ;

            _kernel->evaluateRectangular(configurations.data(), timeRange , {iParticle,iParticle} , range ,  ) ;
        }
    }

    private:

    std::shared_ptr<kernel2B> _kernel;
    int setA;
    int setB;


};