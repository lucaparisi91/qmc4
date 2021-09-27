
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

using range_t = std::array<int,2> ;
#define TOL 1e-7


TEST( greenFunction , caoBerne )
{   
    Real a = 1; Real timeStep=0.1;
    pimc::caoBernePropagator G(timeStep,a);
    std::array<Real,3> x1= {1.3,1.5,0.05};
    std::array<Real,3> x2= {1.6, 1.7,0.1};

    
    Real testValue=G.logEvaluate( x1 , x2  );
    ASSERT_NEAR( testValue,-1.3692611723001526e-06 , 1e-12);

    //x1= {2,0,0};
    //x2= {-2, 0,0};

    testValue=G.logTimeDerivative(x1,x2);

    ASSERT_NEAR(testValue,-0.0001802941377974785, 1e-9);

    testValue=G.logGradientLeft(x1,x2,0);
    
    ASSERT_NEAR(testValue,1.2298012794103977e-05, 1e-9);

    testValue=G.logGradientRight(x1,x2,0);

    ASSERT_NEAR(testValue,9.241786667146941e-06, 1e-9);
    
}

TEST( kernels, pairProductKernel )
{   
    int N=10;
    Real beta = 1;
    int nBeads= 10;
    Real a = 1; Real timeStep=0.1;

    #if DIMENSIONS == 3

    auto G = std::make_shared<pimc::caoBernePropagator>(timeStep,a);

    pimc::pairProductKernel<pimc::caoBernePropagator> pairKernel(G);

    Eigen::Tensor<double,3> data(N,DIMENSIONS,nBeads);
    pimc::geometryPBC_PIMC geo({3000,3000,3000});

    
    pairKernel.setGeometry(geo);
    int t0=4;
    int l=3;
    pairKernel.evaluateTriangular( data , {t0,t0+l} , { N/2 , N - 1} , {N/2, N - 1}  );


    #endif
}

TEST( kernels, primitiveApproximationTwoBodyKernel )
{   
    int N=10;
    Real beta = 1;
    int nBeads= 10;
    Real a = 1; Real timeStep=0.1;

    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );    


    pimc::primitiveApproximationTwoBodyKernel<decltype(V)> pairKernel(std::make_shared<decltype(V)>(V));

    Eigen::Tensor<double,3> data(N,DIMENSIONS,nBeads);
    pimc::geometryPBC_PIMC geo({3000,3000,3000});
    data.setRandom();
    pimc::mask_t mask(N,nBeads);


    pairKernel.setGeometry(geo);
    pairKernel.setTimeStep(timeStep);


    
    
    range_t particleRangeA { 0 , N/2 };
    range_t particleRangeB { N/2+1 , N - 1 };

    range_t timeRange { 4, 7};

    auto sumKernel = pairKernel.evaluateTriangular( data , timeRange , particleRangeA , particleRangeA  );



    Real sumCheck=0;
    for(int t=timeRange[0];t<=timeRange[1];t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== timeRange[1] )) ? 0.5 : 1;
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++)
        {   
            for(int j=0;j<i;j++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);

                }

                sumCheck+=0.5*prefix*r2;
            }
        }
    }

    ASSERT_NEAR(sumCheck*timeStep , sumKernel, TOL);

}


TEST_F( configurationsTest , primitiveApproximation )
{   
    int N=20;
    Real beta = 1;
    int nBeads= 10;
    //int seed = time(NULL);
    int seed = 14;

    SetUp({N/2,N/2},nBeads,1.0);
    SetSeed(seed);
    
    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );    
    

    auto pairKernel = std::make_shared<pimc::primitiveApproximationTwoBodyKernel<decltype(V)> >(std::make_shared<decltype(V)>(V));

    pairKernel->setGeometry(geo);
    pairKernel->setTimeStep(timeStep);

    auto S = std::make_shared<pimc::actionTwoBody>();
    S->setKernel(pairKernel);
    S->setTimeStep(timeStep);
    S->setGeometry(geo);
    S->setSets({0,0});

    SetRandom();

    configurations.fillHeads();

    const auto & data = configurations.dataTensor();

    // test evaluation
    range_t timeRange={0,nBeads-1};
    range_t particleRangeA=configurations.getGroups()[0].range();
    range_t particleRangeB=configurations.getGroups()[1].range();


    Real sumCheck=0;
    for(int t=timeRange[0];t<=timeRange[1]+1;t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== timeRange[1]+1 )) ? 0.5 : 1;
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++)
        {   
            for(int j=particleRangeA[0];j<i;j++)
            {
                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);

                }

                sumCheck+=0.5*prefix*r2;
            }
        }
    }

     auto sum = S->evaluate(configurations,timeRange,particleRangeA);

     ASSERT_NEAR(sum,sumCheck*timeStep,TOL);

    sum = S->evaluateTimeDerivative(configurations,timeRange,particleRangeA);

    ASSERT_NEAR(sum,sumCheck,TOL);



    timeRange=range_t{3,7};
    auto particleRange=range_t{ particleRangeB[0]+2,particleRangeB[0]+7};


    sumCheck=0;
    for(int t=timeRange[0];t<=(timeRange[1]+1);t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== (timeRange[1]+1 ))) ? 0.5 : 1;
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {   
            for(int j=particleRangeA[0];j<i;j++)
            {    
                    Real r2=0;
                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);
                    }

                    sumCheck+=0.5*prefix*r2;
                }
            }
    }

     for(int t=timeRange[0];t<=(timeRange[1]+1);t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== (timeRange[1]+1 ))) ? 0.5 : 1;
        for(int i=particleRange[1]+1;i<=particleRangeA[1];i++)
        {   
            for(int j=particleRange[0];j<=particleRange[1];j++)
            {    
                    Real r2=0;
                    for(int d=0;d<DIMENSIONS;d++)
                    {
                        r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);
                    }

                    sumCheck+=0.5*prefix*r2;
                }
            }
    }



    sum = S->evaluate(configurations,timeRange,particleRange);

    ASSERT_NEAR(sum,sumCheck*timeStep,1e-9);

    //--------- check forces A = B -------------
    timeRange=range_t{0,nBeads-1};

    Eigen::Tensor<Real,3> forces(particleRangeB[1]+1, DIMENSIONS,nBeads);
    Eigen::Tensor<Real,3> forcesCheck(particleRangeB[1]+1, DIMENSIONS,nBeads);

    forces.setConstant(0);
    forcesCheck.setConstant(0);

    for(int t=timeRange[0];t<=timeRange[1];t++){ 
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++) { 
            for(int j=particleRangeA[0];j<i;j++) {
                Real r2=0;
                std::array<Real,3> diff;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    diff[d]=geo.difference(data(i,d,t)-data(j,d,t),d);
                    r2+=diff[d]*diff[d];
                }
                Real r = std::sqrt(r2);
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real tmp = r * diff[d]/r  ;
                    forcesCheck(i,d,t)+=tmp;
                    forcesCheck(j,d,t)-=tmp;
                }
                
                
            }
        }
    }

    S->addGradient(configurations,{0,nBeads-1},particleRangeA,forces);

    for(int t=timeRange[0];t<=timeRange[1]+1;t++){ 
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++) { 
            for(int d=0;d<DIMENSIONS;d++)
            {
                ASSERT_NEAR(forcesCheck(i,d,t),forces(i,d,t),TOL);
            }
                
            }
    }


    //-------------- A different from B ---------------------
    S->setSets({0,1});
    timeRange=range_t{3,7};
    particleRange=range_t{2,7};
    
    sum = S->evaluate(configurations,timeRange,particleRange);

    sumCheck=0;
    for(int t=timeRange[0];t<=timeRange[1]+1;t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== timeRange[1]+1 )) ? 0.5 : 1;
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {   
            for(int j=particleRangeB[0];j<=particleRangeB[1];j++)
            {

                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);
                }

                sumCheck+=0.5*prefix*r2;    
            }
        }
    }

    ASSERT_NEAR(sum,sumCheck*timeStep, 1e-9);

    timeRange=range_t{3,7};
    particleRange=range_t{particleRangeB[0]+2,particleRangeB[0]+7};

    sum = S->evaluate(configurations,timeRange,particleRange);


    sumCheck=0;
    for(int t=timeRange[0];t<=timeRange[1]+1;t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== timeRange[1]+1 )) ? 0.5 : 1;
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {   
            for(int j=particleRangeA[0];j<=particleRangeA[1];j++)
            {

                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);
                }
                sumCheck+=0.5*prefix*r2;    
            }
        }
    }
    ASSERT_NEAR(sum,sumCheck*timeStep,TOL);

    sum = S->evaluateTimeDerivative(configurations,timeRange,particleRange);

    ASSERT_NEAR(sum,sumCheck,TOL);



    //--------- check forces A /= B -------------
    timeRange=range_t{0,nBeads-1};

    forces.setConstant(0);
    forcesCheck.setConstant(0);

    for(int t=timeRange[0];t<=timeRange[1];t++){ 
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++) { 
            for(int j=particleRangeB[0];j<=particleRangeB[1];j++) {
                Real r2=0;
                std::array<Real,3> diff;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    diff[d]=geo.difference(data(i,d,t)-data(j,d,t),d);
                    r2+=diff[d]*diff[d];
                }
                Real r = std::sqrt(r2);
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real tmp = r * diff[d]/r  ;
                    forcesCheck(i,d,t)+=tmp;
                    forcesCheck(j,d,t)-=tmp;
                }
                
            }
        }
    }


    S->addGradient(configurations,timeRange,particleRangeB,forces);

    for(int t=timeRange[0];t<=timeRange[1];t++){ 
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++) { 
            for(int d=0;d<DIMENSIONS;d++)
            {
                ASSERT_NEAR(forcesCheck(i,d,t),forces(i,d,t),TOL);
            }
                
            }
    }

    // check the grand Canonical
    S->setSets({0,0});
    SetGrandCanonicalEnsamble(0);
    int tHead=8;
    int tTail=2;
    int iWorm=3;
    configurations.setHeadTail(iWorm,tHead,tTail);
    configurations.fillHeads();

    timeRange = range_t {0, M-1};
    sum = S->evaluate(configurations,timeRange,particleRangeA);


    sumCheck=0;
    for(int t=timeRange[0];t<=timeRange[1]+1;t++)
    {
        Real prefix= ((t == timeRange[0]) or (t== timeRange[1]+1 )) ? 0.5 : 1;
        for(int i=particleRangeA[0];i<=particleRangeA[1];i++)
        {   
            for(int j=particleRangeA[0];j<i;j++)
            {
                Real mask=prefix;
                if ((i==iWorm) or (j==iWorm))
                {
                    if ((t <= tTail ) or (t>tHead) )
                    {
                        mask=0;
                    }
                    if (( t == tTail + 1) or (t== tHead) )
                    {
                        mask=0.5;
                    }
                }

                Real r2=0;
                for(int d=0;d<DIMENSIONS;d++)
                {
                    r2+=std::pow( geo.difference(data(i,d,t)-data(j,d,t),d),2);

                }

                sumCheck+=0.5*mask*r2;
            }
        }
    }

    ASSERT_NEAR(sum,sumCheck*timeStep,1e-9);


}
