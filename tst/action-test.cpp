#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"
#include "../pimc/hdf5IO.h"
#include "../pimc/toolsPimcTest.h"
#include <filesystem>
#include "../pimc/forces.h"
#include "testConfigurations.h"
#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"

class testAction : public configurationsTest
{
    public:

    Real evaluateHarmonicOnChain( const Eigen::Tensor<Real,3> & data ,  const pimc::geometryPBC_PIMC & geo, const  std::array<int,2> & timeRange, int iChain )
    {
        Real sum=0;
        for(int t=timeRange[0];t<=timeRange[1];t++)
        {
            Real prefactor = ( ( t == timeRange[0]) or (t == timeRange[1]) ) ? 0.5 : 1;

            Real r2=0;
            for(int d=0;d<getDimensions();d++)
            {
                r2+=data(iChain,d,t) * data(iChain,d,t);
            }
            sum+=prefactor*r2*0.5;
        }
        return sum;
    }



};

TEST_F(testAction,oneBodyGrandCanonical )
{
    int N=100;
    int M=50;
    
    SetUp(N,M,1.0);
    
    configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations.dataTensor().setRandom();

    int iOpen=N-1;
    int newHead=M/2;

    int iChainHead=configurations.pushChain(0);
    ASSERT_EQ(iChainHead,N);

    configurations.setHead(iOpen,M);
    configurations.setHead(iChainHead,newHead);
    configurations.setTail(iChainHead,-1);
    configurations.join(iOpen,iChainHead); 
    configurations.fillHeads();

    TESTCHAIN(configurations,iChainHead,newHead,-1);
    
    pimc::geometryPBC_PIMC geo(300,300,300);

     auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) { return 0.5*r*r ;} ,
         [](Real r) {return r  ;}
         );


    const auto & mask = configurations.getTags();

    const auto & data = configurations.dataTensor();
    
    

    
     Real sum=0;
    for(int t=0;t<=M;t++)
    {
        Real prefactor = ((t ==0) or (t == M)) ? 0.5 : 1;
        for(int i=0;i<N;i++)
        {
            sum+=prefactor*V( TRUNCATE_D(data(i,0,t), data(i,1,t),data(i,2,t) ) );
        }
    }



     for(int t=0;t<=newHead;t++)
    {
        Real prefactor = (t ==0 or t == newHead) ? 0.5 : 1;
        int i=N;
        sum+=prefactor*V( TRUNCATE_D(data(i,0,t), data(i,1,t),data(i,2,t) ) );
    }


    auto sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    auto sumAction=sOneBody->evaluate(configurations);

    ASSERT_NEAR(sumAction,sum*timeStep,1e-5);

}   

TEST_F(testAction,oneBodyGrandCanonicalTail )
{
    int N=2;
    int M=10;
    SetUp(N,M,1.0);
    SetSeed(10);
    SetRandom();
    SetUpNonInteractingHarmonicAction();
    SetGrandCanonicalEnsamble(0 );
    
    int tTail=-1;
    int tHead=4;

    configurations.setHeadTail(0,M,tTail);
    configurations.setHeadTail(1,tHead,-1);
    configurations.join(0,1);

    configurations.fillHeads();

    auto & sPot = S.getPotentialAction();

    Real sum = sPot.evaluate(configurations,{tTail+1,M-1},0);

    auto sumCheck = evaluateHarmonicOnChain(configurations.dataTensor(),geo,{tTail+1,M},0);
    ASSERT_NEAR(sum,sumCheck*timeStep,TOL);

    sumCheck=evaluateHarmonicOnChain(configurations.dataTensor(),geo,{0,tHead},0);
    sum = sPot.evaluate(configurations,{0,tHead-1},0);

    ASSERT_NEAR(sum,sumCheck*timeStep,TOL);

    sumCheck=evaluateHarmonicOnChain( configurations.dataTensor(),geo,{tTail+1,M-2},0);

    sum = sPot.evaluate(configurations,{0,M-3},0);

    ASSERT_NEAR( sum, sumCheck*timeStep, TOL);
    
    sum = sPot.evaluate(configurations);

    sumCheck=evaluateHarmonicOnChain(configurations.dataTensor(),geo,{0,tHead},1);
    sumCheck+=evaluateHarmonicOnChain( configurations.dataTensor(),geo,{tTail+1,M},0);

    ASSERT_NEAR(sum,sumCheck*timeStep,TOL);

    sum = sPot.evaluate(configurations,{3,2},0);

    ASSERT_NEAR(sum,0,TOL);

    tTail=4;
    tHead=0;
    configurations.setHeadTail(0,M,tTail);
    configurations.setHeadTail(1,tHead,-1);
    configurations.join(0,1);
    
    SetRandom();
    configurations.fillHeads();

    sum = sPot.evaluate(configurations,{M-1,M-1},0);
    sumCheck=evaluateHarmonicOnChain(configurations.dataTensor(),geo,{M-1,M},0);

    ASSERT_NEAR(sum,sumCheck*timeStep,TOL);
    sum = sPot.evaluate(configurations,{0,0},0);
    ASSERT_NEAR(sum,0,TOL);


}


TEST_F(testAction,twoBody_grandCanonical)
{
    int N1=10;
    int N2=20;

    int M=100;
    Real Beta = 1;
    int buffer=100;

    int nChains = N1 + N2 + 2*buffer;


    Real timeStep=Beta/M;
    srand(11);

    SetUp({N1,N2},M,Beta);

    SetGrandCanonicalEnsamble(0);

    configurations.dataTensor().setRandom();
    configurations.fillHeads();


    auto & data = configurations.dataTensor();
    Real V0=1;
    Real alpha=1;

    auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );
    

   
    //auto sV=pimc::potentialActionTwoBody<decltype(V)>(timeStep,nChains,M,V ,geo,0,0);

    auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);

    pimc::actionTwoBody sV;
    sV.setSets({0,0});
    sV.setKernel(kernel);


    Real r2=0;
    Real sum=0;


    const auto & groupA=configurations.getGroups()[0];
    const auto & groupB=configurations.getGroups()[1];



    for(int t=0;t<=M;t++)
    {
        Real prefactor= ( (t ==0) or (t==M) ) ? 0.5 : 1;

        for(int i=groupA.iStart;i<=groupA.iEnd;i++)
            for(int j=groupA.iStart;j<i;j++)
            {
                Real r2=0;
                for(int d=0;d<getDimensions();d++)
                {
                    r2+= std::pow( geo.difference( data(i,d,t) - data(j,d,t) ,d)    , 2);
                }
                sum+=prefactor*V0*exp(-alpha*r2);
            }
    }


    auto sum2=sV.evaluate(configurations,{0,M-1},{groupA.iStart,groupB.iEnd});

    ASSERT_NEAR(sum*timeStep,sum2,1e-7);


    int iChainHead=0; int tHead=80 *M/100.;
    int iChainTail=1; int tTail=30 *M/100.;


    configurations.setHead(iChainHead, tHead);
    configurations.setTail(iChainTail, tTail);

    
    for(int t=tHead;t<=M;t++)
    {
        Real prefactor= ( (t ==tHead) or (t==M) ) ? 0.5 : 1;
        int i=iChainHead;
        for(int j=groupA.iStart;j<=groupA.iEnd;j++)
            {
                if (j!=i)
                {                    
                    Real r2=0;
                    for(int d=0;d<getDimensions();d++)
                    {
                        r2+= std::pow( geo.difference( data(i,d,t) - data(j,d,t) ,d)    , 2);
                    }
                    sum-=prefactor*V0*exp(-alpha*r2);
                }
            }
    }

    for(int t=0;t<=(tTail + 1);t++)
    {
        Real prefactor= ( (t == (tTail + 1) ) or (t==0) ) ? 0.5 : 1;
            int i=iChainTail;
            for(int j=groupA.iStart;j<=groupA.iEnd;j++)
            {
                if (j != i)
                {
                    Real r2=0;
                for(int d=0;d<getDimensions();d++)
                {
                    r2+= std::pow( geo.difference( data(i,d,t) - data(j,d,t) ,d)    , 2);
                }
                sum-=prefactor*V0*exp(-alpha*r2);
                }
                
            }
    }

    sum2=sV.evaluate(configurations,{0,M-1},{groupA.iStart,groupA.iEnd});

    ASSERT_NEAR(sum*timeStep, sum2,1e-6);

    //auto sVAB=pimc::potentialActionTwoBody<decltype(V)>(timeStep,nChains,M,V ,geo,0,1);


    auto kernelAB = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernelAB->setTimeStep(timeStep);
    kernelAB->setGeometry(geo);

    pimc::actionTwoBody sVAB;
    sVAB.setSets({0,1});
    sVAB.setKernel(kernelAB);

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA,groupB});



    configurations2.setChemicalPotential({0.,0});
    configurations2.dataTensor().setRandom();
    configurations2.fillHeads();

    const auto  & data2 = configurations2.dataTensor();

    sum=0;
    for(int t=0;t<=M;t++)
    {
        Real prefactor= ( (t ==0) or (t==M) ) ? 0.5 : 1;

        for(int i=groupA.iStart;i<=groupA.iEnd;i++)
            for(int j=groupB.iStart;j<=groupB.iEnd;j++)
            {
                Real r2=0;
                for(int d=0;d<getDimensions();d++)
                {
                    r2+= std::pow( geo.difference( data2(i,d,t) - data2(j,d,t) ,d)    , 2);
                }
                sum+=prefactor*V0*exp(-alpha*r2);
            }
    }

    sum2= sVAB.evaluate(configurations2,{0,M-1},{0,N1-1});

    ASSERT_NEAR(sum*timeStep,sum2,1e-7);

}


TEST_F(testAction,twoBody)
{
    int N=100;
    int M=10;
    Real Beta = 1;

    Real timeStep=Beta/M;
    srand(11);

    pimc::geometryPBC_PIMC geo(300,300,300);

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    configurations.dataTensor().setRandom();
    configurations.fillHeads();


    auto & data = configurations.dataTensor();

    // Test on a rectangular interaction potential


    Real Rc=0.1;
    Real V0=2.;

    auto V = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return (r*r)<= Rc*Rc ? V0 :  0 ;} ,
         [](Real r) {return 0  ;}
         );



    //auto sV=pimc::potentialActionTwoBody<decltype(V)>(timeStep,N,M,V ,geo,0,0);


    auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);

    pimc::actionTwoBody sV;
    sV.setSets({0,0});
    sV.setKernel(kernel);



    int t0=0;
    int t1=M-1;
    int iChain=0;

    Real count=0;
    for(int t=t0;t<=t1+1;t++)
    {
         Real prefactor = 1;

        if ( (t==t0) or (t==t1+1) )
        {
            prefactor=0.5;
        }

            for(int j=0;j<N;j++)
            {
                
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(iChain,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (dis<=Rc*Rc and j!= iChain)
                {
                    count+=prefactor;
                }
            }
    }


    auto currentAction = sV.evaluate(configurations,{t0,t1},iChain);
    ASSERT_NEAR(currentAction,count*V0*timeStep,1e-5);




    count=0;
    for(int t=0;t<=M;t++)
    {
        Real prefactor = 1;

        if ( (t==0) or (t==M) )
        {
            prefactor=0.5;
        }

        for(int i=0;i<N;i++)
            for(int j=0;j<i;j++)
            {
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(i,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (dis<=Rc*Rc)
                {
                    count+=prefactor;
                }
            }
    }

    auto totalAction = sV.evaluate(configurations);



    ASSERT_NEAR(totalAction,count*V0*timeStep,1e-5);

    Eigen::Tensor<Real,3> gradientBuffer(N,getDimensions(), M );
    Eigen::Tensor<Real,3> gradientBufferTest(N,getDimensions(), M );


    gradientBuffer.setConstant(0);
    gradientBufferTest.setConstant(0);
    

    sV.addGradient(configurations,{0,M-1},{0,N-1},gradientBuffer);

    Eigen::Tensor<Real,0> sumSquares = (gradientBuffer*gradientBuffer).sum();

    ASSERT_NEAR(sumSquares(0),0,1e-5);

    V0=1.;
    // test on a gaussian interaction potential
    Real alpha = 2.;


     auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );
    

   
    //auto sV2=pimc::potentialActionTwoBody<decltype(V2)>(timeStep,N,M,V2 ,geo,0,0);

    auto kernel2 = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V2) >  >(std::make_shared<decltype(V2)>(V2));

    kernel2->setTimeStep(timeStep);
    kernel2->setGeometry(geo);

    pimc::actionTwoBody sV2;
    sV2.setSets({0,0});
    sV2.setKernel(kernel2);



    Real currentGaussV=0;
    for(int t=t0;t<=t1;t++)
    {
        Real prefactor = 1;
        if ( (t==0) or (t==(M-1)) )
        {
            prefactor=0.5;
        }

        for(int j=0;j<N;j++)
            {
                Real dis=0;
                for(int d=0;d<getDimensions();d++)
                {
                    Real tmp=geo.difference(data(iChain,d,t)-data(j,d,t) ,d);
                    dis+=tmp*tmp;
                }

                if (j != iChain)
                {
                    currentGaussV+=prefactor*V0*exp(-alpha*dis);
                }
                
            }
    }



    currentAction = sV2.evaluate(configurations,{t0,t1-1},iChain);
    ASSERT_NEAR(currentAction,currentGaussV*timeStep,1e-5);

    currentGaussV=0;
    for(int t=0;t<=M;t++)
    {
         Real prefactor = 1;
        if ( (t==0) or (t==M) )
        {
            prefactor=0.5;
        }
        for(int i=0;i<N;i++)
        for(int j=0;j<i;j++)
            {
                Real dis2=0;
                std::array<Real,3> diff;

                for(int d=0;d<getDimensions();d++)
                {
                    diff[d]=geo.difference(data(i,d,t)-data(j,d,t) ,d);
                    dis2+=diff[d]*diff[d];
                }

                currentGaussV+=prefactor*V0*exp(-alpha*dis2);

            }
    }

        for(int t=0;t<=M-1;t++)
    {
        
        for(int i=0;i<N;i++)
        for(int j=0;j<i;j++)
            {
                Real dis2=0;
                std::array<Real,3> diff;

                for(int d=0;d<getDimensions();d++)
                {
                    diff[d]=geo.difference(data(i,d,t)-data(j,d,t) ,d);
                    dis2+=diff[d]*diff[d];
                }

                for(int d=0;d<getDimensions();d++)
                {
                    gradientBufferTest(i,d,t)+=-2*alpha*diff[d]*V0*exp(-alpha*dis2)*timeStep;
                    gradientBufferTest(j,d,t)-=-2*alpha*diff[d]*V0*exp(-alpha*dis2)*timeStep;
                }

            }
    }

    totalAction = sV2.evaluate(configurations);
    ASSERT_NEAR(totalAction,currentGaussV*timeStep,1e-5);

    sV2.addGradient(configurations,{0,M-1},{0,N-1},gradientBuffer);

    Eigen::Tensor<Real,0> sumSquaresTest = (gradientBufferTest*gradientBufferTest).sum();

    sumSquares = (gradientBuffer*gradientBuffer).sum();

    ASSERT_NEAR(sumSquares(0),sumSquaresTest(0),1e-5);

    
}


#include "../pimc/primitivePropagator.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/kernelConstructor.h"


TEST_F(testAction, pairApproximation_primitive)
{
    int N=2;
    int M=10;
    Real Beta = 1;

    Real timeStep=Beta/M;

    SetSeed(135);
    pimc::geometryPBC_PIMC geo(10,10,10);

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    configurations.dataTensor().setRandom();
    configurations.fillHeads();


    auto & data = configurations.dataTensor();
    Real alpha=1;
    Real V0=1;

    auto V = pimc::makeIsotropicPotentialFunctor(
         [V0, alpha](Real r) {return V0* exp(-alpha*r*r) ;} ,
         [V0 , alpha ](Real r) {return -2*V0*alpha*r*exp(-alpha*r*r) ;} 
         );    

    auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);
    auto sV2B=std::make_shared<pimc::actionTwoBody>();
    sV2B->setSets({0,0});
    sV2B->setKernel(kernel);
    sV2B->setTimeStep(timeStep);
    sV2B->setGeometry(geo);

    using propagator_t=pimc::primitivePropagator<decltype(V)>;


    auto G= std::make_shared< propagator_t >(timeStep,V);

    pimc::pairProductKernelConstructor<propagator_t> pairKernelCreator;

    pairKernelCreator.setTimeStep(timeStep);
    pairKernelCreator.setNBeads(M);
    pairKernelCreator.setNMaxParticles(N);
    pairKernelCreator.setGeometry(geo);
    pairKernelCreator.setGreenFunction(G);


    auto pairKernel = pairKernelCreator.create();


    auto SPP = std::make_shared<pimc::actionTwoBody>();
    SPP->setTimeStep(timeStep);
    SPP->setGeometry(geo);
    SPP->setKernel(pairKernel);
    SPP->setSets({0,0});

    auto sum = sV2B->evaluate(configurations);
    auto sumPP = SPP->evaluate(configurations);

    ASSERT_NEAR(sum,sumPP,1e-6);

    Eigen::Tensor<Real,3> forces(N,DIMENSIONS,M+1);
    Eigen::Tensor<Real,3> forces2(N,DIMENSIONS,M+1);

    forces.setConstant(0);
    forces2.setConstant(0);

    sV2B->addGradient(configurations,{0,M-1},{0,N-1},forces);
    SPP->addGradient(configurations,{0,M-1},{0,N-1},forces2);
    
    for ( int i=0;i<N;i++)
    {
        for(int d=0;d< DIMENSIONS ;d++)
        {
            forces(i,d,0)+=forces(i,d,M);
            forces2(i,d,0)+=forces2(i,d,M);
        }
    }

    for(int t=0;t<M;t++)
    {
        for ( int i=0;i<N;i++)
        {
            for(int d=0;d< DIMENSIONS ;d++)
            {
                ASSERT_NEAR( forces(i,d,t) , forces2(i,d,t) ,1e-6) ;
            
            }
        }
    }
}