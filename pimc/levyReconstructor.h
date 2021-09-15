#ifndef LEVYRECONSTRUCTOR_H
#define LEVYRECONSTRUCTOR_H

#include "pimcConfigurations.h"
#include "traits.h"
#include "toolsPimc.h"
#include "action.h"
namespace pimc
{

    enum chainBoundary {  free=0 , fixed=1 };

    class levyReconstructor
        {
            public : 

            levyReconstructor( int maxReconstructionLength);

            void apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,
            const action & S,randomGenerator_t & randG);
            
           

            void setReconstructorBoundaries(chainBoundary lower,chainBoundary upper );

            Real probabilityInitialPositionGaussian(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const;


            Real probabilityInitialPositionUniform(const geometry_t & geo) const;

            Real probabilityInitialPosition(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const;


            void setSigma(Real sigma){sigma2=sigma*sigma;}

            void setTimePBC(bool doTimePBC){applyTimePBC=doTimePBC;}
            
            void setGaussianParticleSampling()
            {
                uniformDistributionFirstParticle=false;
            }

            void setUniformParticleSampling()
            {
                uniformDistributionFirstParticle=true;
            }
            


            private:
            
             void sampleFirstParticle(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG);



            void sampleFirstParticleGaussian(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG);



            void sampleFirstParticleUniform(const geometry_t & geo , std::array<Real,getDimensions()> & x1 ,randomGenerator_t & randG);





             void applyFixedBoundaries(configurations_t & configurations, std::array<int,2> timeRange,int iChain ,
            const action & S,randomGenerator_t & randG);


            void applyFixedFreeBoundaries(configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S, const std::array<int,2> & currentChainBoundary, randomGenerator_t & randG );



            std::array<Real,getDimensions()> mean;
            std::normal_distribution<Real> gauss;
            Eigen::Tensor< Real , 2 > buffer;
            int _maxReconstructionLength;

            std::array<int,2> chainBoundaries;
            bool applyTimePBC;
            configurationsSampler confsSampler;

            bool uniformDistributionFirstParticle;
            std::uniform_real_distribution<float> uniformRealNumber;
            Real sigma2=0;




        };



};

#endif