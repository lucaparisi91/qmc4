#include "levyReconstructor.h"

namespace pimc
{


    levyReconstructor::levyReconstructor( int maxReconstructionLength) :  gauss(0,1),buffer(maxReconstructionLength*2,getDimensions()) , applyTimePBC(false),uniformDistributionFirstParticle(true) , chainBoundaries{ chainBoundary::fixed,chainBoundary::fixed } ,sigma2(1)
    {

    }

    void levyReconstructor::apply (configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S,randomGenerator_t & randG)
    {


        if (
            (chainBoundaries[0] == chainBoundary::fixed) and
        (chainBoundaries[1] == chainBoundary::fixed) 
            )
            {
                applyFixedBoundaries(configurations,timeRange,iChain,S,randG);
            }

        else if ( (chainBoundaries[0] == chainBoundary::free) xor 
                (chainBoundaries[1] == chainBoundary::free)
           )
           {
                
            applyFixedFreeBoundaries(configurations,timeRange,iChain,S, chainBoundaries,randG);
           }
        else if ( (chainBoundaries[0] == chainBoundary::free) and 
                (chainBoundaries[1] == chainBoundary::free)
           )
           {
            std::array<Real , getDimensions()> firstParticleX;
            auto & data = configurations.dataTensor();

            sampleFirstParticle(S.getGeometry() , firstParticleX, randG );

            for(int d=0;d<getDimensions() ;d++)
            {
                data(iChain,d,timeRange[0])=firstParticleX[d];
           }
            
               applyFixedFreeBoundaries(configurations,timeRange,iChain,S,{chainBoundary::fixed,chainBoundary::free},randG);
           }
        
    }


    Real levyReconstructor::probabilityInitialPositionGaussian(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
            {
                Real sum=0;

                for(int d=0;d<getDimensions();d++)
                {
                    sum+= -0.5*log(2*M_PI*sigma2) -0.5 * x[d]*x[d] / sigma2;    
                }
                return sum;
            }
    
    Real levyReconstructor::probabilityInitialPositionUniform(const geometry_t & geo) const
            {
                Real V=1;
                for(int d=0;d<getDimensions();d++)
                {
                    V*=geo.getLBox(d)  ;  
                }
                return -log(V);
            }
    
   Real levyReconstructor::probabilityInitialPosition(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
    {
        if (not uniformDistributionFirstParticle)
        {
            return probabilityInitialPositionGaussian(geo,x);
        }
        else
        {
            return probabilityInitialPositionUniform(geo);
        }
    }

            

    
    void levyReconstructor::sampleFirstParticle(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG)
            {
                if (uniformDistributionFirstParticle)
                {
                    sampleFirstParticleUniform(geo,x1,randG);
                }
                else
                {
                    sampleFirstParticleGaussian(geo,x1,randG);
                }

             }

    void levyReconstructor::sampleFirstParticleGaussian(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG)
    {
        for(int d=0;d<getDimensions();d++)
        {
            x1[d]=gauss(randG)*std::sqrt(sigma2);
        }
    }



    void levyReconstructor::sampleFirstParticleUniform(const geometry_t & geo , std::array<Real,getDimensions()> & x1 ,randomGenerator_t & randG)
    {
        for(int d=0;d<getDimensions();d++)
        {
            x1[d]=(uniformRealNumber(randG)-0.5)*geo.getLBox(d);
        }

    }




     void levyReconstructor::applyFixedFreeBoundaries(configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S, const std::array<int,2> & currentChainBoundary, randomGenerator_t & randG)
     {

         std::array <Real,getDimensions()> startPosition,sampledPosition;
         auto timeStep = S.getTimeStep();
         int l = timeRange[1] - timeRange[0] ;

         auto & data = configurations.dataTensor();

            Real mass = configurations.getGroupByChain(iChain).mass;

            int tFixedBoundary=-1,tFreeBoundary=-1;

            // figure out which boundary is fixed and which is free
            if (chainBoundaries[0]== chainBoundary::fixed)
            {
                tFixedBoundary=timeRange[0];
                tFreeBoundary=timeRange[1];
            }
            else
            {
                tFixedBoundary=timeRange[1];
                tFreeBoundary=timeRange[0];
            }

            // sample the particle position of the free boundary as a gaussian centered on the opposite boundary

            for (int d=0;d<getDimensions();d++)
            {
                startPosition[d]=data(iChain,d,tFixedBoundary);
            }

            confsSampler.sampleFreeParticlePosition(sampledPosition,startPosition,timeStep*l,randG,mass);

            for (int d=0;d<getDimensions();d++)
            {
                data(iChain,d,tFreeBoundary)=sampledPosition[d];
            }

            // sample internal beads using a fixed boundary condition
            applyFixedBoundaries(configurations,timeRange,iChain,S,randG);

     }

    void levyReconstructor::applyFixedBoundaries(configurations_t & configurations, std::array<int,2> timeRange,int iChain ,const action & S,randomGenerator_t & randG)
    {
        /* Does not take into account time periodic boundary conditions
        Reconstruct between timeRange[1] -1 and timeRange[0] + 1
        timeRange[0] and timeRange[1] are used as boundary conditions
        */

        Real timeStep=S.getTimeStep();
        const auto & geo = S.getGeometry();

        int l = timeRange[1] - timeRange[0];

        auto & data = configurations.dataTensor() ;


        const auto & currentChain = configurations.getChain(iChain);
        int iChainNext=currentChain.next;
        

        //int lastBeadChain= timeRange[1] > configurations.nBeads() ? iChainNext : iChain ;
        //int lastBeadTime= timeRange[1] > configurations.nBeads() ? timeRanges[1][1] : timeRange[1] ;


 /*         if (lastBeadChain == -1 )
        {
            throw invalidState("Crossing a head in reconstruction.");
        }

        for (int d=0;d<getDimensions();d++)
            {
                data(iChain,d,timeRange[1])=data(lastBeadChain,d,lastBeadTime);
            }
  */

        std::array<Real,getDimensions()> xBottom;
        std::array<Real,getDimensions()> xTop;
        
        // performs the actual copy
        for (int t=0;t<l-1;t++)
        {

            for (int d=0;d<getDimensions();d++)
                {

                Real eta = gauss(randG);

                //xBottom[d]=geo.difference( data(iChain,d,t+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 
                
                //xTop[d]=geo.difference( data(iChain,d,l+timeRange[0]) - data(iChain,d,t+1+timeRange[0]),d) + data(iChain,d,t+1+timeRange[0]) ; 

                xBottom[d]=data(iChain,d,t+timeRange[0]);
                xTop[d]=data(iChain,d,l+timeRange[0]);


                mean[d] = 
                (
                xBottom[d]*(l-t-1) 
                + xTop[d] )
                /( l - t) 
                ;

                Real variance = (l-t-1) * 1. /(l-t) *timeStep;

                data(iChain,d, timeRange[0] + t + 1 ) = mean[d] + eta *sqrt(variance) ;
                }
        }

        if (applyTimePBC)
        {
            configurations.copyData(  {configurations.nBeads() , timeRange[1] } , iChain , 0  , iChainNext);
        
        }
        
    }







    void levyReconstructor::setReconstructorBoundaries( chainBoundary lower,chainBoundary upper )
    {
        chainBoundaries[0]=lower;
        chainBoundaries[1]=upper;        
    }
       
}