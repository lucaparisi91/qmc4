#include "traits.h"
#include "kernels.h"


namespace pimc
{

template<class greenFunction_t>
class pairProductKernel : public kernel2B
{
    public:
    pairProductKernel( std::shared_ptr<const greenFunction_t> greenFunction_) : 
    greenFunction(greenFunction_) {}
    
    Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB) const
    {
        Real value=0;

        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;
        
        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
            for ( int j=rangeB[0];j<i;j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                    deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);           

                }
                sum+=greenFunction->logEvaluate(deltaX,deltaXNext); 
            }     
        }

        return sum;
    }


    Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB, const mask_t & mask) const
    {
        Real value=0;

        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;

        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                for ( int j=rangeB[0];j<i;j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                        deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);           

                    }
                    if(mask(i,t)*mask(j,t) == 1 )
                    {
                        sum+=greenFunction->logEvaluate(deltaX,deltaXNext) ;
                    }
                     
                }     
        }

        return sum;
    }


    Real evaluateTimeDerivativeTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB) const
    {
        Real value=0;

        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;
        
        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
            for ( int j=rangeB[0];j<i;j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                    deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);           

                }
                sum+=greenFunction->logTimeDerivative(deltaX,deltaXNext); 
            }     
        }

        return sum;
    }
    
    
    Real evaluateRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB) const
    {
        Real value=0;
        
        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;

        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for ( int i=rangeA[0];i<=rangeA[1];i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                    deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);           

                }
                sum+=greenFunction->logEvaluate(deltaX,deltaXNext); 
                
            }     
        }

        return sum;
    }


    Real evaluateTimeDerivativeRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB) const
    {
        Real value=0;
        
        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;

        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for ( int i=rangeA[0];i<=rangeA[1];i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                    deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);           

                }
                sum+=greenFunction->logTimeDerivative(deltaX,deltaXNext); 
            }     
        }

        return sum;
    }

    virtual void addForceRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const
    {
        Real value=0;
        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;


        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for ( int i=rangeA[0];i<=rangeA[1];i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);

                       

                    }

                    for(int d=0;d<getDimensions();d++)
                    {
                        Real left= 
                        greenFunction->logGradientLeft(deltaX,deltaXNext,d);

                        Real right= greenFunction->logGradientRight(deltaX,deltaXNext,d) ;

                        forces(i,d,t)+= left;
                        forces(i,d,t+1)+= right;
                        forces(j,d,t)-= left;
                        forces(j,d,t+1)-= right;

                    }
                    
                }
        }

    }

    virtual void addForceTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) const
    {
       

        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;

        for (int t=timeRange[0];t<=timeRange[1];t++)
        {
            for ( int i=rangeA[0];i<=rangeA[1];i++)
                for ( int j=rangeB[0];j<i;j++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);
                    }


                    for(int d=0;d<getDimensions();d++)
                    {
                        Real left= 
                        greenFunction->logGradientLeft(deltaX,deltaXNext,d);

                        Real right= greenFunction->logGradientRight(deltaX,deltaXNext,d) ;

                        forces(i,d,t)+= left;
                        forces(j,d,t)+= -left;
                        forces(i,d,t+1)+= right;
                        forces(j,d,t+1)+= -right;

                    }
                    
                }
        }

    }

Real evaluateRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, const mask_t & mask) const
    {
        Real value=0;

        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;
        
        Real sum=0;
        for (int t=timeRange[0];t<timeRange[1];t++)
        {
            for ( int i=rangeA[0];i<=rangeA[1];i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    deltaX[d]=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                    deltaXNext[d]=geometry().difference( tn(i,d,t+1) - tn(j,d,t+1) ,d);
                }

                Real ij_mask= mask(i,t) * mask(j,t) ;

                if (ij_mask == 1)
                {
                    sum+=ij_mask*greenFunction->logEvaluate(deltaX,deltaXNext);
                }
                 
            }     
        }

        return sum;
    }


    
    private:


    std::shared_ptr<const greenFunction_t> greenFunction;

};


}