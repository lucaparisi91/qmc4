#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/propagators.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/toolsPimcTest.h"
#include "testConfigurations.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/actionTwoBody.h"
#include "mpi.h"
#include "pairProductKernelOptimized.h"
#include "geometryPMC.h"

Real minimumDistance( const pimc::configurations_t & confs, const std::array<int,2> timeRange, const  std::array<int,2> particleRange, const pimc::geometryPBC_PIMC  & geo )
{
    Real rMin=1e+29;
    const auto & tn = confs.dataTensor();

    for( int t=timeRange[0] ; t<=timeRange[1] ; t++)
    {
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {   
            for(int j=particleRange[0];j<i;j++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diff=geo.difference( tn(i,d,t) -tn(j,d,t) ,d );
                    r2+=diff*diff;
                }
                rMin=std::min(rMin,std::sqrt(r2) );

            }
            
        }
    }

    return rMin;
}

TEST_F( configurationsTest , caoBerne_timing )
{   
    int N=20;
    Real beta=1;
    int nBeads=100;
    Real a=0.001;
    int seed=5678;

    SetUp({N},nBeads,beta);
    

    configurations.fillHeads();


     auto G = std::make_shared<pimc::caoBernePropagator>(beta/nBeads,a);

    pimc::geometry_PBCStorage kernelGeo({1,1,1});

    pimc::geometryPBC_PIMC geoUnrestricted({1,1,1});

    pimc::pairProductKernel<pimc::caoBernePropagator> pairKernel(G);
    pairKernel.setGeometry(geoUnrestricted);

    SetSeed(seed);

    SetRandom();
         while ( minimumDistance(configurations,{0,nBeads},{0,N-1},geo) < a )
        {
            SetRandom();
        }

    auto & data = configurations.dataTensor();
    for(int t=0;t<nBeads+1;t++)
    {
        for(int i=0;i<N;i++)
        {
            for(int d=0;d<DIMENSIONS;d++)
            {
                data(i,d,t)=kernelGeo.pbc(data(i,d,t),d);
            }
            
        }
    }


    int nTrials=1000;
    Real deltat=0;
    Real sum0=0;
    for(int i=0;i<nTrials;i++)
    {
        
        Real t0= MPI_Wtime();
        sum0+=pairKernel.evaluateTriangular( configurations.dataTensor(), {5,84 } , {6,6} , {0,N-1} );
        Real t1= MPI_Wtime();
        deltat+=t1-t0;
    }
    std::cout << std::scientific <<  deltat/nTrials << std::endl;


    pimc::pairProductKernelBufferedDistances<pimc::caoBernePropagator> pairKernelBuffered(G);


    pairKernelBuffered.setGeometry(kernelGeo);

    auto buffer = std::make_shared<Eigen::Tensor<Real,3> >(N*N,DIMENSIONS,nBeads+1);
    buffer->setConstant(0.);
    pairKernelBuffered.setDistanceBuffer(buffer);

    SetSeed(seed);

    deltat=0;
    Real sum1=0;
    for(int i=0;i<nTrials;i++)
    {
        Real t0= MPI_Wtime();
        sum1+=pairKernelBuffered.evaluateTriangular( configurations.dataTensor(), {5,84 } , {6,6} , {0,N-1} );
        Real t1= MPI_Wtime();
        deltat+=t1-t0;
    }


    ASSERT_NEAR(sum0,sum1,1e-9);
    
    std::cout << std::scientific <<  deltat/nTrials << std::endl;


    auto GIsoTropic = std::make_shared<pimc::caoBernePropagator>(beta/nBeads,a);

    pimc::pairProductKernelBufferedDistancesIsotropic<pimc::caoBernePropagator> pairKernelBufferedIsotropic(GIsoTropic);

    pairKernelBufferedIsotropic.setGeometry(kernelGeo);
    auto radiusBuffer = std::make_shared<Eigen::Tensor<Real,2> >(N*N,nBeads+1);
    radiusBuffer->setConstant(0);
    buffer->setConstant(0.);


    pairKernelBufferedIsotropic.setDistanceBuffer(buffer);
    pairKernelBufferedIsotropic.setRadiusBuffer(radiusBuffer);

    
    deltat=0;
    Real sum2=0;
    for(int i=0;i<nTrials;i++)
    {
        
        
        Real t0= MPI_Wtime();
        sum2+=pairKernelBuffered.evaluateTriangular( configurations.dataTensor(), {5,84 } , {6,6} , {0,N-1} );
        Real t1= MPI_Wtime();
        deltat+=t1-t0;
    }
    
    ASSERT_NEAR(sum1,sum2,1e-9);

    
    
    std::cout << std::scientific <<  deltat/nTrials << std::endl;

    std::cout << sum2 << std::endl;


}



TEST_F( configurationsTest , caoBerneGrandCanonicalTiming )
{   
    int N=10;
    Real beta = 1;
    int nBeads= 80;
    //int seed = time(NULL);
    int seed = 14;
    SetUp({N},nBeads,1.0,{2,2,2});
    SetSeed(seed);
    SetGrandCanonicalEnsamble(0);
    SetRandom();

    Real radius=0.01;

     auto G = std::make_shared<pimc::caoBernePropagator>(beta/nBeads,radius);

     
    
     auto kernel = std::make_shared<pimc::pairProductKernel<pimc::caoBernePropagator> >(G);

     kernel->setGeometry(geo);
     kernel->setTimeStep(beta/nBeads);
     
    auto S2B= std::make_shared<pimc::actionTwoBody>();
    S2B->setKernel(kernel);
    S2B->setTimeStep(beta/M);
    S2B->setSets( {0,0});
    S2B->setGeometry(geo);
    S2B->setMinimumDistance(radius);
    

    auto S2B_skipped= std::make_shared<pimc::actionTwoBody>();
    auto policy=std::make_shared<pimc::twoBodySkippingEvaluationPolicy>();
    S2B_skipped->setTwoBodyEvaluationPolicy(policy);
    S2B_skipped->setKernel(kernel);
    S2B_skipped->setTimeStep(beta/M);
    S2B_skipped->setSets( {0,0});
    S2B_skipped->setGeometry(geo);
    S2B_skipped->setMinimumDistance(radius);
    



    while (not S2B->checkConstraints(configurations,{0,   nBeads-1},configurations.getGroups()[0].range() ) )
    {
        SetRandom();
    } ;

    int tTail=4, tHead=30;

    configurations.setHeadTail(3,tHead,tTail);
    //configurations.setHeadTail(8,tHead,-1);
    //configurations.join(3,8);
    configurations.fillHeads();

    Real sum0=0;
    Real dt=0;
    int nTrials=1000;
    int iChain=6;
     for(int i=0;i<nTrials;i++)
    {
        Real t0= MPI_Wtime();
        sum0+=S2B->evaluate(configurations,{10,30} , {iChain,iChain} );
        Real t1= MPI_Wtime();
        dt+=t1-t0;
    }
    std::cout << dt << std::endl;
 
    Real sum1=0;
    dt=0;
    for(int i=0;i<nTrials;i++)
    {
        Real t0= MPI_Wtime();
        sum1+=S2B_skipped->evaluate(configurations,{10,30} , {iChain,iChain} );
        Real t1= MPI_Wtime();
        dt+=t1-t0;
    }
    std::cout << dt << std::endl;

    ASSERT_NEAR(sum0,sum1,1e-9);
}


TEST_F( configurationsTest , caoBerne_accelerator_timing )
{   
    int N=10;
    Real beta=1;
    int nBeads=10;
    Real a=0.01;
    

    int seed=5678;
    std::array<Real,DIMENSIONS> lBox {1,1,1};

    SetUp({N},nBeads,beta,lBox);
    configurations.fillHeads();

     auto G = std::make_shared<pimc::caoBernePropagator>(beta/nBeads,a);


    pimc::geometry_PBCStorage kernelGeo({1,1,1});

  
    pimc::geometryPBC_PIMC geoUnrestricted({lBox[0],lBox[1],lBox[2]});

    SetSeed(seed);

    SetRandom();
         while ( minimumDistance(configurations,{0,nBeads},{0,N-1},geo) < a )
        {
            SetRandom();
        }
    

    auto & data = configurations.dataTensor();
    for(int t=0;t<nBeads+1;t++)
    {
        for(int i=0;i<N;i++)
        {
            for(int d=0;d<DIMENSIONS;d++)
            {
                data(i,d,t)=kernelGeo.pbc(data(i,d,t),d);
            }
            
        }
    }

    std::array<int,2> timeRange {0,1};
    std::array<int,2> particleRange {6,6};


  
    pimc::pairProductKernelBufferedDistances<pimc::caoBernePropagator> pairKernelBuffered(G);
    pairKernelBuffered.setGeometry(kernelGeo);



    auto buffer = std::make_shared<Eigen::Tensor<Real,3> >(N*N,DIMENSIONS,nBeads+1);
    buffer->setConstant(0.);
    pairKernelBuffered.setDistanceBuffer(buffer);


    SetSeed(seed);

    
    Real sum1=pairKernelBuffered.evaluateTriangular( configurations.dataTensor(), timeRange , particleRange , {0,N-1} );

    std::array< int, DIMENSIONS> nCells {2,2,2};


    pimc::cellNeighbourList acc;
    acc.setGeometry(geo);
    acc.setNBeads(M);
    acc.setNCells(nCells);

    acc.build();
    
    acc.add(configurations.dataTensor(),{0,M}, {0,N-1} );

    Real sum2=pairKernelBuffered.evaluateTriangular( configurations.dataTensor(), acc, timeRange , particleRange , {0,N-1} );

    ASSERT_NEAR(sum1,sum2,1e-9);
    

}