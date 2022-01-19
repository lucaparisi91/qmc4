#ifndef PAIR_PRODUCT_KERNEL
#define PAIR_PRODUCT_KERNEL


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

        
        Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB) 
        {
            
            Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            

            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                       
                            
                            sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            }



            return sum;
            
        }

        
        Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB, const mask_t & mask) 
        {
            Real value=0;

            Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                            if(mask(i,t) * mask(j,t) == 1  )
                            {
                                sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                            }
                        k++;
                    }     
            }



            return sum;
        }


        Real evaluateTimeDerivativeTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB) 
        {

                       Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                       
                            
                            sum+=greenFunction->logTimeDerivative( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            }



            return sum;

        }
        
        
        Real evaluateRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB) 
        {
           Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                       
                            
                            sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            }



            return sum;
            
        }


        Real evaluateTimeDerivativeRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB) 
        {
                       Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                       
                            
                            sum+=greenFunction->logTimeDerivative( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            }
            return sum;
        }

        virtual void addForceRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) 
        {
             auto & buffer = distanceBuffer();


            for (int t=timeRange[0];t<=timeRange[1]+1;t++)
            {
                int k=0;
                for ( int i=rangeA[0];i<=rangeA[1];i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                        
                            for(int d=0;d<getDimensions();d++)
                            {
                                
                                buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                            }
                            k++;

                        

                    }
            }

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for ( int i=rangeA[0];i<=rangeA[1];i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            Real left= 
                            greenFunction->logGradientLeft( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)},d);

                            Real right= greenFunction->logGradientRight( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)},d);
                            
                            forces(i,d,t)+= left;
                            forces(j,d,t)+= -left;
                            forces(i,d,t+1)+= right;
                            forces(j,d,t+1)+= -right;
                        }
                        k++;
                    }
            }


        }

        virtual void addForceTriangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, Eigen::Tensor<Real,3> & forces) 
        {

             auto & buffer = distanceBuffer();


            for (int t=timeRange[0];t<=timeRange[1]+1;t++)
            {
                int k=0;
                for ( int i=rangeA[0];i<=rangeA[1];i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        
                            for(int d=0;d<getDimensions();d++)
                            {
                                
                                buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                            }
                            k++;

                        

                    }
            }

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for ( int i=rangeA[0];i<=rangeA[1];i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            Real left= 
                            greenFunction->logGradientLeft( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)},d);

                            Real right= greenFunction->logGradientRight( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)},d);
                            
                            forces(i,d,t)+= left;
                            forces(j,d,t)+= -left;
                            forces(i,d,t+1)+= right;
                            forces(j,d,t+1)+= -right;
                        }
                        k++;
                    }
            }
        }

    Real evaluateRectangular(const Eigen::Tensor<Real,3> & tn, const  std::array<int,2> & timeRange, const std::array<int,2> & rangeA, const std::array<int,2> & rangeB, const mask_t & mask)
        {

            Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                        for(int d=0;d<getDimensions();d++)
                        {
                            
                            buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                        }
                        k++;

                    }     
            }

            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<=rangeB[1];j++)
                    {
                       
                            auto ij_mask= mask(i,t) * mask(j,t) ;

                            if (ij_mask == 1)
                            {
                            sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;

                            }
                        k++;
                    }     
            }

            return sum;

        }

        auto & distanceBuffer() {return *_distances; }
        void  setDistanceBuffer(std::shared_ptr<Eigen::Tensor<Real,3> > buffer) {_distances=buffer; }



        
        private:


        std::shared_ptr<const greenFunction_t> greenFunction;
        std::shared_ptr<Eigen::Tensor<Real,3> > _distances;


    };


}

#endif