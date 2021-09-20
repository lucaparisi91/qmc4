
#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/propagators.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/toolsPimcTest.h"
#include "testConfigurations.h"
#include "../pimc/potentialKernel.h"

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

    pairKernel.setGeometry(geo);
    int t0=4;
    int l=3;

    
    pairKernel.evaluateTriangular( data , {t0,t0+l} , { N/2 , N - 1} , {N/2, N - 1}  );

}

