#ifndef KERNELS_PIMC_H
#define KERNELS_PIMC_H

#include <memory>
#include "pimcConfigurations.h"
#include "geometryPMC.h"

namespace pimc
{

using mask_t = maskTensor;
using configurations_t = pimcConfigurations;


class kernel2B
{
    public:

    kernel2B(){}

    const auto &  geometry () const {return _geometry; } ; 
    
   
    void setGeometry ( const geometryPBC_PIMC & new_geometry_) { _geometry=new_geometry_ ;};


    virtual Real evaluateRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {
        if (configurations.getEnsamble() == ensamble_t::canonical)
        {
            return evaluateRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            return evaluateRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB,configurations.getTags());
        }
        
    }


    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const = 0;


    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn, 
        const  std::array<int,2> & timeRange,
        const std::array<int,2> & rangeA,
        const std::array<int,2> & rangeB,
        const mask_t & mask) const = 0;


    virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const = 0;




    private:
   
    geometryPBC_PIMC _geometry;
   
};

};

#endif