#ifndef GEOMETRYPMC
#define GEOMETRYPMC


#include "geometry.h"
#include <unsupported/Eigen/CXX11/Tensor>

using Real = double;

namespace pimc
{


class geometryPBC_PIMC : public geometryPBC
{
    public:
    using geometryPBC::geometryPBC;

    void  updateSpringDifferences(Eigen::Tensor<Real,3> & differences,  const  Eigen::Tensor<Real, 3> & tn , std::array< int, 2 > timeRange, std::array<int,2>  particleRange    );     // update spring distances stored in an array ( time_slices, dimensions, number of particles)
    

    void updateEqualTimeDifferences( Eigen::Tensor<Real,3>  & differences, const  Eigen::Tensor<Real, 3> & tn, std::array<int,2> timeRange, int i, std::array<int,2> particleRange );



};


struct geometry_PBCStorage
{
    geometry_PBCStorage(){}
    geometry_PBCStorage(std::array<Real,DIMENSIONS> lBox)  {
        setLbox(lBox);
    }

    Real pbc( Real x, int d) const { if (x >= _lBoxInverseHalf[d] ) { x-=_lBox[d]; } else if( x<= -_lBoxInverseHalf[d]) {x+=_lBox[d];} ; return x; }

    Real difference(Real x,int d) const  {return pbc(x,d);}

    void setLbox( std::array<Real,DIMENSIONS> lBox_) {
        _lBox=lBox_;
        for(int d=0;d<DIMENSIONS;d++)
        {
            _lBoxInverseHalf[d]=0.5*1./(_lBox[d]);
        }
    }


    private:

    std::array<Real,DIMENSIONS> _lBox;
    std::array<Real,DIMENSIONS> _lBoxInverseHalf;

};


};

#endif