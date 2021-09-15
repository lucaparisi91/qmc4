#include "twoBodyPotential.h"
#include "utils.h"
#include <cassert>


using Real = double;
using namespace particleKernels;


template<int dims>
twoBodyPotential<dims>::twoBodyPotential(std::array<int,2> rangeA_,std::array<int,2> rangeB_, std::array<Real,dims> lBox_  ,
std::array<int, 3> dimensions) :
_dimensions(dimensions)
{
    setPeriodic(lBox_);
    setDimensions(dimensions);

    setRange(rangeA_,rangeB_);
    
};


template<int dims>
void twoBodyPotential<dims>::postUpdateRange()
{
    int iStartA=rangeA[0];
    int iStartB=rangeB[0];

    int iEndA=rangeA[1];
    int iEndB=rangeB[1];
 
    // check that the two ranges do not intersect
    assert(iStartA<=iEndA);
    assert(iStartB<=iEndB);
    
    if ((iStartA == iStartB) and (iEndA == iEndB) )
    {
        isTriangular=true;
    }
    else
    {
        isTriangular=false;

        int lowerSet= iStartA < iStartB ? 0 : 1;

        if (lowerSet == 0)
        {
            assert(iEndA<iStartB);
            
        }
        else
        {
            assert(iEndB<iStartA);
        }
    }

}


template<int dims>
bool twoBodyPotential<dims>::containedInSetA(int i1,int i2) const
{
    return (i1 >= rangeA[0]) and (i2 <=rangeA[1]);
}

template<int dims>
bool twoBodyPotential<dims>::containedInSetB(int i1,int i2) const
{
    return (i1 >= rangeB[0]) and (i2 <=rangeB[1]);
};


template<int dims>
void twoBodyPotential<dims>::setPeriodic(std::array<Real,dims> lBox   )
{
    periodic=true;
    for(int d=0;d<dimensions;d++)
    {
        _lBox[d]=lBox[d];
    }

};


template<int dims>
void twoBodyPotential<dims>::setParticleDimensions(std::array<int,3> dimensions   )
{
    _dimensions=dimensions;
};

template<int dims>
void twoBodyPotential<dims>::setForceDimensions(std::array<int,3> dimensions   )
{
    _forceDimensions=dimensions;
};


template<int dims>
void twoBodyPotential<dims>::setDimensions(std::array<int,3> dimensions   )
{
    setParticleDimensions(dimensions);
    setForceDimensions(dimensions); 
};

template<int dims>
void twoBodyPotential<dims>::setRange( const std::array<int,2> & rangeA_,const std::array<int,2> & rangeB_   )
{
    rangeA=rangeA_;
    rangeB=rangeB_;

    postUpdateRange();
};

template class twoBodyPotential<3>;
template class twoBodyPotential<2>;
template class twoBodyPotential<1>;


