#ifndef PAIRPRODUCT_MESH_KERNEL_H
#define PAIRPRODUCT_MESH_KERNEL_H

#include "particleContainer/particleDistancesList.h"

namespace pimc
{
    struct meshKernel
    {
        virtual Real evaluate(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, int iParticle )=0;

        virtual Real evaluate(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, const range_t & particleRange )=0;

        virtual Real evaluateTimeDerivative(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, const range_t & particleRange )=0;


        virtual void addForce(const Eigen::Tensor<Real,3> & tn, const linkedCellParticles & particles, const  std::array<int,2> & timeRange, const std::array<int,2> & particleRange, Eigen::Tensor<Real,3> & forces)=0;


    };


    template<class greenFunction_t>
    class pairProductMeshKernel : public meshKernel 
    {
        public:

        pairProductMeshKernel( std::shared_ptr<greenFunction_t> greenFunction_ , size_t N, size_t M) :
        _greenFunction(greenFunction_),
        disBuffer(N,M)
        {
            
        }

        void setGeometry(const geometry_t & geo) {_geo=geo;}


        Real evaluate(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, int iParticle )
        {
            disBuffer.buildDistanceList(tn,particles,_geo,{timeRange[0],timeRange[1] },iParticle);


            const auto & buffer = disBuffer.getBufferDistances();
            const auto & bufferOff =  disBuffer.getBufferDistancesOffset() ;
            const auto & nDistances = disBuffer.getNDistances();
            

            Real sum=0;

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                for(int k=0;k<nDistances[t];k++)
                    {
                        sum+=_greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{bufferOff(k,0,t),bufferOff(k,1,t),bufferOff(k,2,t)}) ;
                    }   
                //std::cout << "sum2: "<< t << " " << sum << std::endl;

            }

            
            return sum;

        }


         Real evaluate(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, const range_t & particleRange )
        {
            disBuffer.buildDistanceList(tn,particles,_geo,{timeRange[0],timeRange[1] },particleRange);

            const auto & buffer = disBuffer.getBufferDistances();
            const auto & bufferOff =  disBuffer.getBufferDistancesOffset() ;
            const auto & nDistances = disBuffer.getNDistances();

            Real sum=0;

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                for(int k=0;k<nDistances[t];k++)
                    {
                        sum+=_greenFunction->logEvaluate( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{bufferOff(k,0,t),bufferOff(k,1,t),bufferOff(k,2,t)}) ;
                    }   
                //std::cout << "sum2: "<< t << " " << sum << std::endl;

            }

            
            return sum;

        }

        Real evaluateTimeDerivative(const Eigen::Tensor<Real,getDimensions()> & tn, const linkedCellParticles & particles,const std::array<int,2> & timeRange, const range_t & particleRange )
        {
            
           disBuffer.buildDistanceList(tn,particles,_geo,{timeRange[0],timeRange[1] },particleRange);

            const auto & buffer = disBuffer.getBufferDistances();
            const auto & bufferOff =  disBuffer.getBufferDistancesOffset() ;
            const auto & nDistances = disBuffer.getNDistances();

            Real sum=0;

            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                for(int k=0;k<nDistances[t];k++)
                    {
                        sum+=_greenFunction->logTimeDerivative( {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)},{bufferOff(k,0,t),bufferOff(k,1,t),bufferOff(k,2,t)}) ;
                    }   
                //std::cout << "sum2: "<< t << " " << sum << std::endl;

            }

            
            return sum;

        }

        virtual void addForce(const Eigen::Tensor<Real,3> & tn, const linkedCellParticles & particles, const  std::array<int,2> & timeRange, const std::array<int,2> & particleRange, Eigen::Tensor<Real,3> & forces) 
        {
            disBuffer.buildDistanceList(tn,particles,_geo,{timeRange[0],timeRange[1] },particleRange);

            const auto & buffer = disBuffer.getBufferDistances();
            const auto & bufferOff =  disBuffer.getBufferDistancesOffset() ;
            const auto & nDistances = disBuffer.getNDistances();

            const auto & pairs = disBuffer.getPairList();


            for (int t=timeRange[0];t<=timeRange[1];t++)
            {
                for(int k=0;k<nDistances[t];k++)
                    {
                        auto i = pairs(k,t,0);
                        auto j = pairs(k,t,1);                        

                        std::array<Real,3> delta {buffer(k,0,t),buffer(k,1,t),buffer(k,2,t)};
                        std::array<Real,3> deltaNext {bufferOff(k,0,t),bufferOff(k,1,t),bufferOff(k,2,t)};


                        for(int d=0;d<getDimensions();d++)
                        {
                            Real left= 
                            _greenFunction->logGradientLeft(  { delta[0],delta[1],delta[2]},{deltaNext[0], deltaNext[1],deltaNext[2]},d);

                            Real right= _greenFunction->logGradientRight(  { delta[0],delta[1],delta[2]},{deltaNext[0], deltaNext[1],deltaNext[2]},d);
                            
                            forces(i,d,t)+= left;
                            forces(j,d,t)+= -left;
                            forces(i,d,t+1)+= right;
                            forces(j,d,t+1)+= -right;
                        }


                    }   

            }

            
        }



        const auto & getDistanceBuffer() { return disBuffer;}


        private:

        std::shared_ptr<greenFunction_t> _greenFunction;
        twoBodyDistancesBuffer disBuffer;
        geometry_t _geo;

    };

}

#endif