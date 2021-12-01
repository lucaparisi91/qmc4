#include "actionTwoBody.h"
namespace pimc
{


Real twoBodyEvaluationPolicy::evaluateTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) 
{
    if (configurations.getEnsamble() == ensamble_t::canonical)
        {
            return _kernel->evaluateTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            return _kernel->evaluateTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB,configurations.getTags() );
        }
}


Real twoBodyEvaluationPolicy::evaluateRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) 
{
    if (configurations.getEnsamble() == ensamble_t::canonical)
        {
            return _kernel->evaluateRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            return _kernel->evaluateRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB,configurations.getTags() );
        }
}


Real twoBodyEvaluationPolicy::evaluateTimeDerivativeRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) 
    {
        if (configurations.getEnsamble() == ensamble_t::canonical or (not configurations.isOpen()) )
        {
            return _kernel->evaluateTimeDerivativeRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented for a grandCanonical open ensamble");
        }

    }


Real twoBodyEvaluationPolicy::evaluateTimeDerivativeTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) 
    {
        if (configurations.getEnsamble() == ensamble_t::canonical or (not configurations.isOpen()) )
        {
            return _kernel->evaluateTimeDerivativeTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented for a grandCanonical open ensamble");
        }

    }



 Real twoBodySkippingEvaluationPolicy::evaluateTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) 
    {
        Real sum=0;

        if (configurations.isOpen(setA()) )
        {
            std::array<int,2> wormChains;
            int nWorms=0;

            wormChains[0]=configurations.getGroups()[setA()].heads[0];
            wormChains[1]=configurations.getGroups()[setA()].tails[0];

            if (wormChains[0] > wormChains[1] )
            {
                std::swap(wormChains[0],wormChains[1]);
            }
            if (wormChains[0]==wormChains[1])
            {
                nWorms=1;
            }
            else
            {
                nWorms=2;
            }

            int j0=rangeB[0];

            for( int i=0 ;i< nWorms ;i++ )
            {
                int iWorm = wormChains[i];
                auto tHead=configurations.getChain(iWorm).head;
                auto tTail=configurations.getChain(iWorm).tail;
                
                auto rangeBB = intersectRanges(rangeB,{j0,iWorm-1});
                rangeBB=intersectRanges(rangeBB,{rangeB[0],rangeA[0]-1});
                



                sum+=kernel().evaluateRectangular( configurations.dataTensor(), timeRange,rangeA, rangeBB);
                

                if (iWorm < rangeA[0])
                {
                    sum+=kernel().evaluateRectangular( configurations.dataTensor(), {std::max(timeRange[0],tTail+1), std::min(tHead,timeRange[1])},rangeA, {iWorm,iWorm});
                }

                j0=iWorm+1;
            }

            sum+=kernel().evaluateTriangular( configurations.dataTensor(), timeRange,rangeA, {j0,rangeB[1]});


        }
        else
        {
              sum+=kernel().evaluateTriangular( configurations.dataTensor(), timeRange,rangeA, rangeB);
        }

        
    return sum;
    } 

/*
Real twoBodySkippingEvaluationPolicy::evaluateRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {
        Real sum=0;

        if (configurations.isOpen(setA()) )
        {
            std::array<int,2> wormChains;

            wormChains[0]=configurations.getGroups()[setA()].heads[0];
            wormChains[1]=configurations.getGroups()[setA()].tails[0];

            if (wormChains[0] > wormChains[1] )
            {
                std::swap(wormChains[0],wormChains[1]);
            }

            int j0=rangeB[0];

            for(auto iWorm : wormChains )
            {
               
                auto tHead=configurations.getChain(iWorm).head;
                auto tTail=configurations.getChain(iWorm).tail;
                
                auto rangeBB = intersectRanges(rangeB,{j0,iWorm-1});
                sum+=kernel().evaluateRectangular( configurations.dataTensor(), {timeRange[0], std::min(tTail,timeRange[1])},rangeA, rangeBB);
                
                sum+=kernel().evaluateRectangular( configurations.dataTensor(), {std::max(timeRange[0],tHead+1), timeRange[1]},rangeA, rangeBB);

                rangeBB = intersectRanges(rangeB,{j0,iWorm});

                sum+=kernel().evaluateRectangular( configurations.dataTensor(), {std::max(timeRange[0],tTail+1), std::min(tHead,timeRange[1])},rangeA, rangeBB);


                j0=iWorm+1;
            }

        }

        
    return sum;
    }

*/

Real actionTwoBody::evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange_, const range_t & particleRange )
    {
        const auto & groupA = configurations.getGroups()[setA];
        const auto & groupB = configurations.getGroups()[setB];

        if (timeRange_[1] < timeRange_[0] )
        {
            return 0;
        }

         range_t timeRange {timeRange_[0],timeRange_[1]+1};

        Real sum=0;

        if (setA == setB)
        {
            auto range=intersectRanges(particleRange, groupA.range());

            sum+=evaluationPolicy->evaluateTriangular(configurations , timeRange, range , groupA.range() );
            sum+=evaluationPolicy->evaluateRectangular(configurations , timeRange, {range[1]+1,groupA.range()[1]} , range );

        }
        else
        {
            auto rangeA = intersectRanges(  groupA.range(),particleRange );
            auto rangeB = intersectRanges(  groupB.range(),particleRange );

            sum+=evaluationPolicy->evaluateRectangular(configurations , timeRange, rangeA , groupB.range() );
            sum+=evaluationPolicy->evaluateRectangular(configurations , timeRange, rangeB , groupA.range() );


        }

        return sum;
    }


Real actionTwoBody::evaluateTimeDerivative(const configurations_t & configurations,const std::array<int,2> & timeRange_, const range_t & particleRange )
    {
        const auto & groupA = configurations.getGroups()[setA];
        const auto & groupB = configurations.getGroups()[setB];

        if (timeRange_[1] < timeRange_[0] )
        {
            return 0;
        }
        
        range_t timeRange {timeRange_[0],timeRange_[1]+1};

        Real sum=0;

        if (setA == setB)
        {
             auto range=intersectRanges(particleRange, groupA.range());

            sum+=evaluationPolicy->evaluateTimeDerivativeTriangular(configurations , timeRange, range , groupA.range() );
            sum+=evaluationPolicy->evaluateTimeDerivativeRectangular(configurations , timeRange, {range[1]+1,groupA.range()[1]} , range );

        }
        else
        {
             auto rangeA = intersectRanges(  groupA.range(),particleRange );
            auto rangeB = intersectRanges(  groupB.range(),particleRange );

            sum+=evaluationPolicy->evaluateTimeDerivativeRectangular(configurations , timeRange, rangeA , groupB.range() );
            sum+=evaluationPolicy->evaluateTimeDerivativeRectangular(configurations , timeRange, rangeB , groupA.range() );

        }

        return sum;
    }





void actionTwoBody::addGradient(const configurations_t & configurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer)
{
    assert( ! configurations.isOpen()  );

    auto const & groupA = configurations.getGroups()[setA];
    auto const & groupB = configurations.getGroups()[setB];

    
    if (setA == setB)
    {
        auto range=intersectRanges(particleRange, groupA.range());

        _kernel->addForceTriangular( configurations.dataTensor(), timeRange, range ,  groupA.range(), gradientBuffer);
        //_kernel->addForceRectangular( configurations.dataTensor(), timeRange, {range[1]+1,groupA.range()[1]} ,  range, gradientBuffer);

    }
    else
    {
        auto rangeA = intersectRanges(  groupA.range(),particleRange );
        auto rangeB = intersectRanges(  groupB.range(),particleRange );

        _kernel->addForceRectangular(configurations.dataTensor(), timeRange , rangeA, groupB.range(),  gradientBuffer);

        _kernel->addForceRectangular(configurations.dataTensor(), timeRange , rangeB, groupA.range(),  gradientBuffer);


    }


}


bool actionTwoBody::checkConstraints(const configurations_t & configurations,const std::array<int,2> & timeRange_,const  std::array<int,2> & particleRange)
{
    if (_distanceMinimumConstraint)    
    {
         if (timeRange_[1] < timeRange_[0] )
        {
            return true;
        }
        

        range_t timeRange {timeRange_[0],timeRange_[1]+1};


        auto const & groupA = configurations.getGroups()[setA];
        auto const & groupB = configurations.getGroups()[setB];
    
        if (setA == setB)
        {
            auto rangeA = intersectRanges(  groupA.range(),particleRange );
            
            bool pass = _kernel->checkMinimumDistanceConstraintTriangular( configurations,timeRange,rangeA , groupA.range(), rMin );

            if (pass)
            {
                pass= _kernel->checkMinimumDistanceConstraintRectangular( configurations,timeRange,rangeA,{ rangeA[1]+1,groupA.range()[1] }, rMin );
            }

            return pass;

        }
        else
        {
             auto rangeA = intersectRanges(  groupA.range(),particleRange );
             auto rangeB = intersectRanges(  groupB.range(),particleRange );

             bool pass = _kernel->checkMinimumDistanceConstraintTriangular( configurations,timeRange,rangeA , groupB.range(), rMin );

            if (pass)
            {
                pass= _kernel->checkMinimumDistanceConstraintRectangular( configurations,timeRange,rangeB,groupA.range(), rMin );
            }

            return pass;

            
        }
    }
    else {

        return true;
    }
}



};

