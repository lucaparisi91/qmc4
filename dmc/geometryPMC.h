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

};

#endif