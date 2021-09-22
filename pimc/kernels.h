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
    const auto &  timeStep () const {return _timeStep; } ;
    void setTimeStep (Real timeStep_)  { _timeStep=timeStep_; } ;
   
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
            return evaluateRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB,configurations.getTags() );
        }

    }

    virtual Real evaluateTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {
        if (configurations.getEnsamble() == ensamble_t::canonical)
        {
            return evaluateTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            return evaluateTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB,configurations.getTags() );
        }

    }


    virtual Real evaluateTimeDerivativeRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {
        if (configurations.getEnsamble() == ensamble_t::canonical or (not configurations.isOpen()) )
        {
            return evaluateTimeDerivativeRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented for a grandCanonical open ensamble");
        }

    }

    virtual Real evaluateTimeDerivativeTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {
        if (configurations.getEnsamble() == ensamble_t::canonical or (not configurations.isOpen()) )
        {
            return evaluateTimeDerivativeTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB);
        }
        else
        {
            throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented for a grandCanonical open ensamble");
        }

    }






    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const =0;


    virtual Real evaluateTimeDerivativeRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const {throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented");}

    virtual Real evaluateTimeDerivativeTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const {throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented");}




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

     virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB,
        const mask_t & mask
    ) const { throw std::runtime_error("evaluateTriangular masked is not implemented"); };

     virtual void addForceRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const = 0 ;

     virtual void addForceTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const = 0;


    private:
   
    geometryPBC_PIMC _geometry;
    Real _timeStep;
};

};

#endif