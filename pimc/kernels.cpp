#include "kernels.h"

namespace pimc
{
    bool kernel2B::checkMinimumDistanceConstraintRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA,const std::array<int,2> & rangeB, Real minDistance)
    {
        Real minDistanceSquared=minDistance*minDistance;

        for (int t=timeRange[0];t<=timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
                {
                    Real r2=0;
                    for(int d=0;d<getDimensions();d++)
                    {
                        Real tmp = geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        r2+=tmp*tmp;
                    }

                    if (r2 <= minDistanceSquared)
                    {
                        return false;
                    }

                }   
        }
        return true;
    } 

    bool kernel2B::checkMinimumDistanceConstraintTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA,const std::array<int,2> & rangeB, Real minDistance)
    {
        Real minDistanceSquared=minDistance*minDistance;
        
        for (int t=timeRange[0];t<=timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                for ( int j=rangeB[0];j<i;j++)
                {
                    Real r2=0;
                    for(int d=0;d<getDimensions();d++)
                    {
                        Real tmp = geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        r2+=tmp*tmp;
                    }

                    if (r2 <= minDistanceSquared)
                    {
                        return false;
                    }

                }   
        }
        return true;
    }

    bool kernel2B::checkMinimumDistanceConstraintTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA,const std::array<int,2> & rangeB, Real minDistance,const mask_t & mask)
    {
        Real minDistanceSquared=minDistance*minDistance;
        
        for (int t=timeRange[0];t<=timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                for ( int j=rangeB[0];j<i;j++)
                {
                    if (mask(i,t)*mask(j,t) == 1 or mask(i,t-1)*mask(j,t-1) == 1 )
                    {

                    
                        Real r2=0;
                        for(int d=0;d<getDimensions();d++)
                        {
                            Real tmp = geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                            r2+=tmp*tmp;
                        }

                        if (r2 <= minDistanceSquared)
                        {
                            return false;
                        }
                    }


                }   
        }
        return true;
    } 


     bool kernel2B::checkMinimumDistanceConstraintRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA,const std::array<int,2> & rangeB, Real minDistance,const mask_t & mask)
    {
        Real minDistanceSquared=minDistance*minDistance;

        for (int t=timeRange[0];t<=timeRange[1];t++)
        {
            for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                for ( int j=rangeB[0];j<=rangeB[1];j++)
                {

                    Real r2=0;
                    if (mask(i,t)*mask(j,t) == 1 or mask(i,t-1)*mask(j,t-1) == 1 )
                    {

                        for(int d=0;d<getDimensions();d++)
                        {
                            Real tmp = geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                            r2+=tmp*tmp;
                        }

                        if (r2 <= minDistanceSquared)
                        {
                            return false;
                        }

                    }

                }   
        }
        return true;
    }
    
        
}



