#include "actionTwoBody.h"
namespace pimc
{

using range_t = std::array<int,2>;


auto intersectRanges(const range_t & rangeA,const range_t & rangeB )
{
    range_t intersection{std::max(rangeA[0],rangeB[0]),std::min(rangeA[1],rangeB[1]) };

    if (intersection[1] < intersection[0])
    {
        return range_t{0,-1};
    }
    else
    {
        return intersection;
    }
}

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

            sum+=_kernel->evaluateTriangular(configurations , timeRange, range , groupA.range() );
            sum+=_kernel->evaluateRectangular(configurations , timeRange, {range[1]+1,groupA.range()[1]} , range );

        }
        else
        {
            auto rangeA = intersectRanges(  groupA.range(),particleRange );
            auto rangeB = intersectRanges(  groupB.range(),particleRange );

            sum+=_kernel->evaluateRectangular(configurations , timeRange, rangeA , groupB.range() );
            sum+=_kernel->evaluateRectangular(configurations , timeRange, {groupA.range()[0], rangeA[0] - 1} , rangeB );
            sum+=_kernel->evaluateRectangular(configurations , timeRange, { rangeA[1] + 1,groupA.range()[1]} , rangeB );

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

            sum+=_kernel->evaluateTimeDerivativeTriangular(configurations , timeRange, range , groupA.range() );
            sum+=_kernel->evaluateTimeDerivativeRectangular(configurations , timeRange, {range[1]+1,groupA.range()[1]} , range );

        }
        else
        {
            auto rangeA = intersectRanges(  groupA.range(),particleRange );
            auto rangeB = intersectRanges(  groupB.range(),particleRange );

            sum+=_kernel->evaluateTimeDerivativeRectangular(configurations , timeRange, rangeA , groupB.range() );
            sum+=_kernel->evaluateTimeDerivativeRectangular(configurations , timeRange, {groupA.range()[0], rangeA[0] - 1} , rangeB );
            sum+=_kernel->evaluateTimeDerivativeRectangular(configurations , timeRange, { rangeA[1] + 1,groupA.range()[1]} , rangeB );
            
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


        _kernel->addForceRectangular(
            configurations.dataTensor()    , timeRange ,
             {groupA.range()[0], rangeA[0] - 1}    
            , rangeB,  gradientBuffer);
        
        _kernel->addForceRectangular(
            configurations.dataTensor()    , timeRange ,
            {rangeA[1]+1, groupA.range()[1] }    
            , rangeB,  gradientBuffer);

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
            
            bool pass = _kernel->checkMinimumDistanceConstraintTriangular( configurations.dataTensor(),timeRange,rangeA , groupA.range(), rMin );

            if (pass)
            {
                pass= _kernel->checkMinimumDistanceConstraintRectangular( configurations.dataTensor(),timeRange,rangeA,{ rangeA[1]+1,groupA.range()[1] }, rMin );
            }

            return pass;

        }
        else
        {
            auto rangeA = intersectRanges(  groupA.range(),particleRange );
            auto rangeB = intersectRanges(  groupB.range(),particleRange );


            bool pass= _kernel->checkMinimumDistanceConstraintRectangular( configurations.dataTensor(),timeRange,rangeA,groupB.range(), _distanceMinimumConstraint );

            pass = pass and _kernel->checkMinimumDistanceConstraintRectangular( configurations.dataTensor(),timeRange,rangeB,groupA.range(), _distanceMinimumConstraint );
            return pass;   
        }
    }
    else {

        return true;
    }
}

};

