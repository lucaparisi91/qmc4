#include <array>
#include "utils.h"


namespace particleKernels
{

template<int dims>
struct twoBodyPotential
{
    using Real=double;
    static const int dimensions = dims;
    
    twoBodyPotential(
        std::array<int,2> rangeA,std::array<int,2> rangeB, // particle index range for set A and B. If A=B particles are consider indistinguishible
        std::array<Real,dims> lBox_ , // length of the box
        std::array<int,3> dimensions // dimensions of the array
         );

    void setPeriodic( std::array<Real,dims> lBox);
    void setRange(const  std::array<int,2> & rangeA_, const std::array<int,2> &  rangeB_);
    

    void setDimensions( std::array<int,3> dimensions);
    void setParticleDimensions( std::array<int,3> dimensions);
    void setForceDimensions( std::array<int,3> dimensions);


    template<class V_t>
    Real operator()(const V_t & V,const Real * positions, int i1, int i2, int t0 , int t1) const;

    template<class V_t>
    Real operator()(const V_t & V,const Real * positions, int i1, int i2, int t0 , int t1,const Real * mask) const;
    

    template<class V_t>
    void addForce(const V_t & V,const Real * positions,Real * forces, int i1 , int i2, int t0 , int t1,Real multiplicationFactor=1) const;


    bool containedInSetA(int i1, int i2) const;
    bool containedInSetB(int i1, int i2) const;



    private:

    void postUpdateRange();


    bool periodic=false;
    std::array<int,2> rangeA;
    std::array<int,2> rangeB;
    std::array<int,3> _dimensions;
    std::array<int,3> _forceDimensions;

    std::array<Real,dimensions> _lBox;
    bool isTriangular=false;
};


#include "twoBodyPotential.hpp"


}


