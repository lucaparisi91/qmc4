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
        for (int t=timeRange[0];t<=timeRange[1];t++)
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


    Real evaluateRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB) const
    {
        Real value=0;
        
        std::array<Real,DIMENSIONS> deltaX;
        std::array<Real,DIMENSIONS> deltaXNext;

        Real sum=0;
        for (int t=timeRange[0];t<=timeRange[1];t++)
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
    

    private:

    std::shared_ptr<const greenFunction_t> greenFunction;

};


}