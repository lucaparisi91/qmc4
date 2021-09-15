#include "gtest/gtest.h"
#include "twoBodyPotential.h"
#include <random>
#include <cmath>
#include "timers.h"
#include "testUtils.h"


using Real = double;
using namespace particleKernels;

TEST(twoBodyPotential, evaluate)
{
    const int N=100;
    const int T=10;

    constexpr int D = 1;
    Real l =1;


    std::array<Real,D> lBox;

    for(int d=0;d<D;d++)
    {
        lBox[d]=l;
    }
    

    double * particles = new double[N*D*T];
    double * forces = new double[N*D*T];
    double * forcesCheck = new double[N*D*T];

    int nTrials = 20;

    setConstant(forces,0,N-1,0,T-1,N,D,T, 0 );
    setConstant(forcesCheck,0,N-1,0,T-1,N,D,T, 0 );



    std::array<int,2> rangeA { 0,N/2};
    std::array<int,2> rangeB { N/2 + 1, N-1};
    std::array<int,2> rangeC { 0 , N-1};
    
    
    int seed=10;
    std::default_random_engine generator(seed);

    generateRandomParticlesInBox(particles,0,N-1,0,T-1,N ,D, T, generator, lBox.data() );
    

    twoBodyPotential<D> potAB(rangeA,rangeB,lBox,{N,D,T});
    twoBodyPotential<D> potAA(rangeA,rangeA,lBox,{N,D,T});
    twoBodyPotential<D> potBB(rangeB,rangeB,lBox,{N,D,T});
    twoBodyPotential<D> potCC(rangeC,rangeC,lBox,{N,D,T});



    potAA.setRange(rangeA,rangeA);

    gaussianPotential V(1);

    Real sum=0;
    Real sumCheck=0;

    for(int nT=0;nT<nTrials;nT++)
    {
        sum = potAB(V,particles,rangeA[0],rangeA[1],0,T-1);
        sumCheck=potAB(V,particles,rangeA[0],rangeB[1],0,T-1);

        ASSERT_NEAR(sumCheck,sum,1e-6);



    



        sumCheck=0;
        for(int t=0;t<=T-1;t++)
            for(int i=rangeA[0];i<=rangeA[1];i++)
                for(int j=rangeB[0];j<=rangeB[1];j++)
                {
                    Real r2=0;
                    
                    for (int d=0;d<D;d++)
                    {
                        Real diffd=particles[ i + N*d + t*D*N] - particles[ j + N*d + t*D*N];
                        diffd=utils::differencePBC(diffd,lBox[d],1./lBox[d]);
                        r2+=diffd*diffd;
                        

                    }
                    sumCheck+=V(std::sqrt(r2));

                }

        ASSERT_NEAR(sum,sumCheck,1e-6);


    }
   
    sum = potAA(V,particles,rangeA[0],rangeA[1],0,T-1);



    sumCheck=0;
    for(int t=0;t<=T-1;t++)
        for(int i=rangeA[0];i<=rangeA[1];i++)
            for(int j=rangeA[0];j<i;j++)
            {
                Real r2=0;
                for (int d=0;d<D;d++)
                {
                    Real diffd=particles[ i + N*d + t*D*N] - particles[ j + N*d + t*D*N];
                    diffd=utils::differencePBC(diffd,lBox[d],1./lBox[d]);
                    r2+=diffd*diffd;
                }
                sumCheck+=V(std::sqrt(r2));
            }
    ASSERT_NEAR(sum,sumCheck,1e-6);

    // test on multiple partial updates



    std::uniform_int_distribution<int> disIndexA(rangeA[0],rangeA[1]);
    std::uniform_int_distribution<int> disIndexB(rangeB[0],rangeB[1]);
    std::uniform_int_distribution<int> disIndexAB(rangeA[0],rangeB[1]);


    auto sumABOld = potAB(V,particles,rangeA[0],rangeA[1],0,T-1);
    auto sumAAOld = potAA(V,particles,rangeA[0],rangeA[1],0,T-1);

    int t0=0;
    int t1=T-1;

    for (int nT=0;nT<nTrials;nT++)
    {
            //std::cout << nT << std::flush << std::endl;
          
            int iStart= disIndexAB(generator);
            int iEnd= disIndexAB(generator);

            auto deltaSumAB=-potAB(V,particles,iStart,iEnd,0,T-1);
            auto deltaSumAA=-potAA(V,particles,iStart,iEnd,0,T-1);
            
            generateRandomParticlesInBox(particles,iStart,iEnd,t0,t1,N ,D,  T, generator, lBox.data() );
            deltaSumAB+=potAB(V,particles,iStart,iEnd,0,T-1);
            deltaSumAA+=potAA(V,particles,iStart,iEnd,0,T-1);

            auto sumABNew = sumABOld +deltaSumAB;
            auto sumAANew = sumAAOld +deltaSumAA;

            auto sumABNewCheck =  potAB(V,particles,rangeA[0],rangeB[1],0,T-1);
            auto sumAANewCheck =  potAA(V,particles,rangeA[0],rangeA[1],0,T-1);

            ASSERT_NEAR(sumABNewCheck,sumABNew,1e-5);
            ASSERT_NEAR(sumAANewCheck,sumAANew,1e-5);

            
            sumAAOld=sumAANewCheck;
            sumABOld=sumABNewCheck;

    }

    {
    // triangular two body forces and potential of a set C should be the same when split into two different sets A + B
    auto sumAA = potAA(V,particles,rangeA[0],rangeB[1],0,T-1);
    auto sumAB = potAB(V,particles,rangeA[0],rangeB[1],0,T-1);
    auto sumBB = potBB(V,particles,rangeA[0],rangeB[1],0,T-1);
    auto sumCC = potCC(V,particles,rangeA[0],rangeB[1],0,T-1);

    sum=sumAA + sumAB + sumBB;

    ASSERT_NEAR(sumCC,sum,1e-6);

    potAA.addForce(V,particles,forces,rangeA[0],rangeA[1],0,T-1);
    potBB.addForce(V,particles,forces,rangeB[0],rangeB[1],0,T-1);
    potAB.addForce(V,particles,forces,rangeA[0],rangeB[1],0,T-1);
    potCC.addForce(V,particles,forcesCheck,rangeA[0],rangeB[1],0,T-1);


    for(int i=0;i<N*D*T;i++)
    {
        ASSERT_NEAR(forcesCheck[i],forces[i],1e-6);
    }



    }




}