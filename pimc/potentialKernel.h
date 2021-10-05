#include "kernels.h"
/*
Describe the primitive approximation two body kernel
*/

namespace pimc
{
template<class V_t>
class primitiveApproximationTwoBodyKernel : public kernel2B
{
    public:

    
    primitiveApproximationTwoBodyKernel(const  std::shared_ptr<V_t> & V_) : V(V_) {};

    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const
    {

    Real sum2b=0;
//#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=timeRange[0];t<=timeRange[1];t++)
        for (int iParticle=rangeA[0];iParticle<=rangeA[1];iParticle++) 
        {
            for (int jParticle=rangeB[0];jParticle<=rangeB[1];jParticle++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diffd=geometry().difference( tn(iParticle,d,t) - tn(jParticle,d,t) ,d);
                    r2+= diffd * diffd;
                }

                Real prefactor= ( (t== timeRange[0]) or (t == timeRange[1]) ) ? 0.5 : 1 ;
                sum2b+=(*V)(std::sqrt(r2))*prefactor;
            }
        }
        return sum2b*timeStep();
    }

    virtual Real evaluateTimeDerivativeRectangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  const
    {
        return evaluateRectangular(tn,timeRange,rangeA,rangeB)/timeStep();
    }


    virtual Real evaluateTimeDerivativeTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  const
    {
        return evaluateTriangular(tn,timeRange,rangeA,rangeB)/timeStep();
    }

    virtual Real evaluateRectangular(
        const Eigen::Tensor<Real,3> & tn, 
        const  std::array<int,2> & timeRange,
        const std::array<int,2> & rangeA,
        const std::array<int,2> & rangeB,
        const mask_t & mask) const
    {

    Real sum2b=0;
//#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=timeRange[0];t<=timeRange[1];t++)
        for (int iParticle=rangeA[0];iParticle<=rangeA[1];iParticle++) 
        {
            for (int jParticle=rangeB[0];jParticle<=rangeB[1];jParticle++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diffd=geometry().difference( tn(iParticle,d,t) - tn(jParticle,d,t) ,d);
                    r2+= diffd * diffd;
                }

                Real ij_mask=  
                        0.5*( mask(iParticle,t) * mask(jParticle,t) + mask(iParticle,t-1) * mask(jParticle,t-1) ) ;
                Real prefactor= ( (t== timeRange[0]) or (t == timeRange[1]) ) ? 0.5 : 1 ;

                ij_mask = ((t == timeRange[0]) or (t==timeRange[1])) ? std::min(ij_mask,0.5) : ij_mask ;

                sum2b+=(*V)(std::sqrt(r2))*ij_mask;
            }
        }
        return sum2b*timeStep();

    }


   virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB,
        const mask_t & mask
    ) const 
    {

    Real sum2b=0;
//#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=timeRange[0];t<=timeRange[1];t++)
        for (int iParticle=rangeA[0];iParticle<=rangeA[1];iParticle++) 
        {
            for (int jParticle=rangeB[0];jParticle< iParticle;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diffd=geometry().difference( tn(iParticle,d,t) - tn(jParticle,d,t) ,d);
                    r2+= diffd * diffd;
                }

                Real ij_mask=  
                        0.5*( mask(iParticle,t) * mask(jParticle,t) + mask(iParticle,t-1) * mask(jParticle,t-1) ) ;
                

                ij_mask = ((t == timeRange[0]) or (t==timeRange[1])) ? std::min(ij_mask,0.5) : ij_mask ;

                sum2b+=(*V)(std::sqrt(r2))*ij_mask;
            }
        }
        return sum2b*timeStep();

    }

    virtual Real evaluateTriangular(
        const Eigen::Tensor<Real,3> & tn
    ,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const 
    {

    Real sum2b=0;
//#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=timeRange[0];t<=timeRange[1];t++)
        for (int iParticle=rangeA[0];iParticle<=rangeA[1];iParticle++)
        {
            for (int jParticle=rangeB[0];jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diffd=geometry().difference( tn(iParticle,d,t) - tn(jParticle,d,t), d );
                    r2+= diffd * diffd;
                }

                Real prefactor= ( (t== timeRange[0]) or (t == timeRange[1]) ) ? 0.5 : 1 ;
                sum2b+=(*V)(std::sqrt(r2))*prefactor;
            }
        }
        return sum2b*timeStep();
    }

   virtual void addForceRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const
     {
//#pragma omp parallel for schedule(static) collapse(3)
        for (int t=timeRange[0];t<=timeRange[1];t++)
            for (int i=rangeA[0];i<=rangeA[1];i++)
            {
                for (int j=rangeB[0];j<=rangeB[1];j++)
                {
                    Real r2=0;
                    std::array<Real,3> diff;

                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        diff[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        r2+= diff[d] * diff[d];
                    }

                    auto r = std::sqrt(r2);
                    auto rInverse = 1./r;
                    auto dVdr = V->radialDerivative(r);                
                    
                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        Real tmp=dVdr*diff[d]*rInverse*timeStep();
                        forces(i,d,t)+=tmp;
                        forces(j,d,t)-=tmp;
                    }
                }
            
            }
    }

     virtual void addForceTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const
     {
//#pragma omp parallel for schedule(static) collapse(3)
        for (int t=timeRange[0];t<=timeRange[1];t++)
            for (int i=rangeA[0];i<=rangeA[1];i++)
            {
                for (int j=rangeB[0];j< i ;j++)
                {
                    Real r2=0;
                    std::array<Real,3> diff;

                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        diff[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        r2+= diff[d] * diff[d];
                    }

                    auto r = std::sqrt(r2);
                    auto rInverse = 1./r;
                    auto dVdr = V->radialDerivative(r);                
                    
                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        Real tmp=dVdr*diff[d]*rInverse*timeStep();
                        forces(i,d,t)+=tmp;
                        forces(j,d,t)-=tmp;
                    }
                }
            
            }
    }

    

    private:
    std::shared_ptr<V_t> V;

};

}