#ifndef PAIR_PRODUCT_KERNEL_OPTIMIZED
#define PAIR_PRODUCT_KERNEL_OPTIMIZED

#include "../pimc/accelerationStructure.h"

namespace pimc
{
    template<class greenFunction_t>
    class pairProductKernelBufferedDistances
    {
        using geometry_t = geometry_PBCStorage;

        public:

        pairProductKernelBufferedDistances( std::shared_ptr<const greenFunction_t> greenFunction_) : 
        greenFunction(greenFunction_) {}
    
    
    const auto & geometry() { return _geometry;}

    void setGeometry ( const geometry_t & new_geometry_) { _geometry=new_geometry_ ;};


    Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB )
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
                       
                            
                            sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            }



            return sum;
        }

   Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, cellNeighbourList & cellList, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB )
        {
            throw std::runtime_error("Not implemented yet");
            Real value=0;
 
            Real sum=0;

            auto & buffer = distanceBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                auto & currentList = cellList.cellNeighBourList(t);

                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                {
                    int iCell=currentList.cellIndex(i);

                    for(int jjCell=0;jjCell<currentList.nCellNeighbours(iCell);jjCell++)
                    {
                        int jCell=currentList.neighbourCell(iCell,jjCell);
                        for ( int jj=0;jj< currentList.nParticles(jCell);jj++)
                        {
                            auto j=currentList.particle(jCell,jj);
                            //std::cout << j << " " << t << std::endl;

                            if (j<i)
                            {
                                int k = j + i*tn.dimensions()[0] ;
                                for(int d=0;d<getDimensions();d++)
                                {
                                    buffer(k,d,t)=geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);

                                }   
                            }
                        }
                    }
                }
            }

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                auto & currentList = cellList.cellNeighBourList(t);

                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                {
                    int iCell=currentList.cellIndex(i);

                    for(int jjCell=0;jjCell<currentList.nCellNeighbours(iCell);jjCell++)
                    {
                        int jCell=currentList.neighbourCell(iCell,jjCell);
                        for ( int jj=0;jj< currentList.nParticles(jCell);jj++)
                        {
                            auto j=currentList.particle(jCell,jj);
                            std::cout << j << " " << t << std::endl;

                            if (j<i)
                            {
                                int k = j + i*tn.dimensions()[0] ;
                                sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                                k++;
                                k++;
                            }

                        }
                    }
                }
            }


            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                const auto & currentList = cellList.cellNeighBourList(t);


                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                {
                    int iCell=currentList.cellIndex(i);

                    for(int jjCell=0;jjCell<currentList.nCellNeighbours(iCell);jjCell++)
                    {
                        int jCell=currentList.neighbourCell(iCell,jjCell);
                        for ( int jj=0;jj< currentList.nParticles(jCell);jj++)
                        {
                            auto j=currentList.particle(jCell,jj);

                            if (j<i)
                            {
                                
                            }
                            

                        }
                    }
                }
            }
            

            
           /*  for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                       
                            
                            sum+=greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{buffer(k,0,t+1),buffer(k,1,t+1),buffer(k,2,t+1)}) ;
                        
                        k++;

                    }     
            } */



            return sum;
        }
    

    auto & distanceBuffer() {return *_distances; }
    void  setDistanceBuffer(std::shared_ptr<Eigen::Tensor<Real,3> > buffer) {_distances=buffer; }
    

    protected:

    std::shared_ptr<Eigen::Tensor<Real,3> > _distances;
    std::shared_ptr<const greenFunction_t> greenFunction;
    geometry_t _geometry;

    };

}


namespace pimc
{
    template<class greenFunction_t>
    class pairProductKernelBufferedDistancesIsotropic :     public pairProductKernelBufferedDistances<greenFunction_t>
    {
        public:

        pairProductKernelBufferedDistancesIsotropic( std::shared_ptr<const greenFunction_t> greenFunction_) : 
        pairProductKernelBufferedDistances<greenFunction_t>::pairProductKernelBufferedDistances(greenFunction_) {}

        Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB )
        {
            Real value=0;

            Real sum=0;

            auto & buffer = this->distanceBuffer();
            auto & radius = radiusBuffer();

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {
                        radius(k,t)=0;
                        for(int d=0;d<getDimensions();d++)
                        {   
                            buffer(k,d,t)=this->geometry().difference( tn(i,d,t) - tn(j,d,t) ,d);
                            radius(k,t)+=buffer(k,d,t)*buffer(k,d,t);
                        }
                        radius(k,t)=std::sqrt(radius(k,t) );

                        k++;

                    }     
            }
            
            for (int t=timeRange[0];t<timeRange[1];t++)
            {
                int k=0;
                for (int i=rangeA[0]; i <= rangeA[1] ; i++)
                    for ( int j=rangeB[0];j<i;j++)
                    {                   
                        Real cosTeta=0;
                        for(int d=0;d<DIMENSIONS;d++)
                        {
                            cosTeta+=buffer(k,d,t)*buffer(k,d,t+1);
                        }
                        cosTeta/=(radius(k,t)*radius(k,t+1));

                        sum+=this->greenFunction->logEvaluate( radius(k,t),radius(k,t+1),cosTeta) ;
                        
                        k++;

                    }     
            }

            return sum;
        }
    

   

    auto & radiusBuffer() {return *_radius; }
    void  setRadiusBuffer(std::shared_ptr<Eigen::Tensor<Real,2> > buffer) {_radius=buffer; }
    protected:
    std::shared_ptr<Eigen::Tensor<Real,2> > _radius;

    };

}

/* 
 template<class greenFunction_t>
    class pairProductKernelCellNeighbourList
    {
        using geometry_t = geometry_PBCStorage;

        public:

        pairProductKernelBufferedDistances( std::shared_ptr<const greenFunction_t> greenFunction_) : 
        greenFunction(greenFunction_) {}
    
    
    const auto & geometry() { return _geometry;}

    void setGeometry ( const geometry_t & new_geometry_) { _geometry=new_geometry_ ;};

    Real evaluateTriangular(const Eigen::Tensor<Real,3> & tn, const std::array<int,2> & timeRange, const std::array<int ,2 > & rangeA , const std::array<int , 2 > & rangeB )
        {

        }
    private:

    std::shared_ptr<const greenFunction_t> greenFunction;
    geometry_t _geometry;


    };
 */
#endif