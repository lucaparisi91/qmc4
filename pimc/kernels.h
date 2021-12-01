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
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  =0;


    virtual Real evaluateTimeDerivativeRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  {throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented");}

    virtual Real evaluateTimeDerivativeTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) {throw std::runtime_error("evaluateTimeDerivativeRectangular not implemented");}




    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn, 
        const  std::array<int,2> & timeRange,
        const std::array<int,2> & rangeA,
        const std::array<int,2> & rangeB,
        const mask_t & mask) = 0;


    virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  = 0;

     virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB,
        const mask_t & mask
    ) { throw std::runtime_error("evaluateTriangular masked is not implemented"); };

     virtual void addForceRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces)  = 0 ;

     virtual void addForceTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces)  = 0;

    virtual bool checkMinimumDistanceConstraintRectangular(const configurations_t & configurations , const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance)
    {
         if (configurations.getEnsamble() == ensamble_t::canonical  )
         {
             return checkMinimumDistanceConstraintRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB,minDistance );
         }
         else
        {
            return checkMinimumDistanceConstraintRectangular(configurations.dataTensor(),timeRange,rangeA,rangeB,minDistance,configurations.getTags() );
        }

    }

    virtual bool checkMinimumDistanceConstraintTriangular(const configurations_t & configurations, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance)
    {
         if (configurations.getEnsamble() == ensamble_t::canonical  )
         {
             return checkMinimumDistanceConstraintTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB,minDistance );
         }
         else
        {
            return checkMinimumDistanceConstraintTriangular(configurations.dataTensor(),timeRange,rangeA,rangeB,minDistance,configurations.getTags() );
        }

    }

    


    virtual bool checkMinimumDistanceConstraintRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance);    

    virtual bool checkMinimumDistanceConstraintTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance);  

    virtual bool checkMinimumDistanceConstraintRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance,const mask_t & mask);    

    virtual bool checkMinimumDistanceConstraintTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Real minDistance, const mask_t & mask);  
    

    
    private:
    
    geometryPBC_PIMC _geometry;
    Real _timeStep;
};

};

#endif