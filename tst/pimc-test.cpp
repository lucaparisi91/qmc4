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

namespace fs = std::filesystem;


TEST(distances,updateSingleParticleDistances)
{
    pimc::geometryPBC_PIMC geo(10,10,10);
    int D = getDimensions();
    int N= 1;

    state_t positions(N,D);
    positions.setRandom();
    positions=positions - 0.5 ;

    auto difference = geo.differencesTwoBody(positions);
    positions(0,0)+=1;
    auto difference2 = geo.differencesTwoBody(positions);

    geo.differencesTwoBody(difference, positions , 0 );

    for (int k=0;k<difference2.rows();k++)
    {
        ASSERT_NEAR( difference2(k,0) , difference(k,0) , 1e-6  ) ;
    }

    int T = 50;

    Eigen::Tensor<Real, 3> data(N,getDimensions(),T+1);
    data.setRandom();


    Eigen::Tensor<Real,3> springDifferences(T,getDimensions() , N );


    geo.updateSpringDifferences(  springDifferences, data, {0,T-1} , {0,N-1} );

    for (int t=0;t<T;t++)
    {
        for (int i=0;i<N;i++)
        {
            int d=0;
            ASSERT_NEAR( springDifferences(t,d,i) , data(i,d,t+1) - data(i,d,t) , 1e-4 );
        }
    }

    Eigen::Tensor<Real,3> potentialDifferences(N ,getDimensions() , T + 1 );

    int jChain = N-1;
    geo.updateEqualTimeDifferences(potentialDifferences,data, {0,T-1}, jChain , {0,N-1});

    for (int t=0;t<T;t++)
    {
        for (int i=0;i<N;i++)
        {
            for(int d=0;d<getDimensions();d++)
            {
            ASSERT_NEAR( potentialDifferences(i,d,t) , data(i,d,t) - data(jChain,d,t) ,1e-4);
            }
        }
    }

    //auto springDistances = geo.springDistances( timeConfigurations );

    //geo.springDistances(springDistances,timeConfigurations,0,2 , 0, 3  );

}

#if DIMENSIONS == 3
TEST(configurations, init)
{
    const int N = 20;
    const int M = 30;

    pimc::particleGroup groupA{ 0 , N-1, N-1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();

    



    pimc::geometryPBC_PIMC geo(10,10,10);

    Real timeStep = 1e-2;



    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads()  , geo);




    Real currentKineticAction = sT.evaluate(configurations);

    Real kineticActionSimple=0;

    for(int n=0;n<N;n++)
    {
        for (int d=0;d<getDimensions();d++)
        {
             ASSERT_NEAR(data(n,d,M), data(n,d,0)  , 1e-5);
        }
        for (int t=0;t<M;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    kineticActionSimple+=std::pow(data(n,d,(t+1)%M ) - data(n,d,t) ,2) /(4*0.5*timeStep) ;   
                }
            }
    }


    ASSERT_NEAR(currentKineticAction,kineticActionSimple,1e-5);

      #if DIMENSIONS == 3

      auto harmonicPotential = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );
    #endif

      #if DIMENSIONS == 1
     auto harmonicPotential = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return x  ;},
         );
    #endif

    {
    Real x=0.1;
    Real y=0.23;
    Real z=0.56;
    
    ASSERT_NEAR(0.5*(x*x + y*y + z*z), harmonicPotential(x,y,z),1e-4);

    }

     pimc::potentialActionOneBody<decltype(harmonicPotential)> pot1(timeStep,harmonicPotential,geo);

     auto v = pot1.evaluate(configurations);
     Real vCheck=0;

     for (int t=1;t<=M-1;t++)
        for(int n=0;n<N;n++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                vCheck+=0.5*data(n,d,t) * data(n,d,t) ;
            }
        }
    
      for(int n=0;n<N;n++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                vCheck+=0.25*data(n,d,0) * data(n,d,0) ;
                vCheck+=0.25*data(n,d,M) * data(n,d,M) ;
            }
        }

        vCheck*=timeStep;
    ASSERT_NEAR(v,vCheck,1e-4);



}

#endif

#if DIMENSIONS == 3


TEST(moves,levy_reconstructor)
{
    int seed = 30;
    const int N = 50;
    const int M = 30;

    
    pimc::geometryPBC_PIMC geo(10,10,10);


    pimc::particleGroup groupA{ 0 , N-1, N ,  1.0};
    randomGenerator_t randG(seed+5);

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    Real timeStep=1e-3;

    std::array<int,2> timeSlice= {10,26};

    pimc::kineticAction sT(timeStep, configurations.nChains(), configurations.nBeads()  , geo);
     

    pimc::levyReconstructor levy(M);
    int iChain = 30;

    configurations.dataTensor().setRandom();
    configurations.fillHeads();

    configurations2=configurations;


    levy.apply(configurations,timeSlice,iChain,sT,randG);

    auto &  data = configurations.dataTensor();
    auto &  data2 = configurations2.dataTensor();

     for (int i=0;i<N;i++)
    {
        for (int t=0; t<=timeSlice[0];t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_NEAR( data2(i,d,t) , data(i,d,t) ,1e-4);
        }

        for (int t=timeSlice[1]; t<M;t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_NEAR( data2(i,d,t) , data(i,d,t) ,1e-4);
        }

        if (i == iChain)
        {
            for (int t=timeSlice[0]+1; t<timeSlice[1];t++)
         {
            for (int d=0;d<getDimensions();d++)
            ASSERT_FALSE( abs( data2(i,d,t) - data(i,d,t) ) < 1e-5);
        }

        };

    }

    // test if gaussian distribution is sampled for a 3 point time slice

    int iSampled=29;
    iChain=3;
    timeSlice={iSampled-1,iSampled+1};
    Real alpha=1./(4*0.5*M_PI*timeStep);

    std::array<Real,getDimensions()> averagePosition {0,0,0};
    averagePosition={0,0,0};

    std::array<Real,getDimensions()> averagePositionSquared {0,0,0};
    std::array<Real,getDimensions()> averagePositionExtremes {0,0,0};
    int nSteps=100000;

    for( int d=0;d<getDimensions();d++)
    {
        averagePositionExtremes[d]=0.5*( data(iChain,d,iSampled-1) + data(iChain,d,(iSampled+1)%M) );
    }

    for (int i=0;i< nSteps ; i++)
    {
        levy.apply(configurations,timeSlice,iChain,sT,randG);


        for(int d=0;d<getDimensions();d++)
        {
            averagePosition[d]+=data(iChain,d,iSampled);
            averagePositionSquared[d]+=data(iChain,d,iSampled) * data(iChain,d,iSampled) ;
        }
    }


    for (int d=0;d<getDimensions();d++)
    {
        averagePosition[d]/=nSteps;
        averagePositionSquared[d]/=nSteps;

        Real sigma= averagePositionSquared[d] - averagePosition[d]*averagePosition[d];
        Real sigmaExpected = timeStep/2.;

        Real meanDeviation=std::abs((averagePosition[d]  -  averagePositionExtremes[d])/averagePositionExtremes[d]);
        Real sigmaDeviation =  std::abs(sigmaExpected - sigma)/sigmaExpected;

        ASSERT_LT(sigmaDeviation,1e-2);
        ASSERT_LT(meanDeviation,1e-3);

    }


}

#endif

TEST(moves,levy)
{
    int seed = 30;
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N, 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});
     auto & data = configurations.dataTensor();
    data.setRandom();
    configurations.fillHeads();

    pimc::geometryPBC_PIMC geo(10,10,10);

    randomGenerator_t randG(seed);

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupA});

    
    Real timeStep=1e-2;
    std::array<int,2> timeSlice= {10,26};

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, N + 1 , M  , geo);
    
    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );    
    #endif

    #if DIMENSIONS == 3
    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );    
    #endif


    std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    pimc::firstOrderAction S(sT,  sV);

    pimc::levyReconstructor levy(M);

    pimc::levyMove mover(20,0);
    int success = 0;

    auto & data2 = configurations2.dataTensor();

    int nSteps=10000;
    for (int i=0;i<nSteps;i++)
    {
        configurations2=configurations;
        bool accept = mover.attemptMove(configurations,S,randG);

        if (accept)
        {
            success+=1;
        }
        else
        {
            for(int t=0;t<M;t++)
                for(int i=0;i<N;i++)
                {
                    for(int d=0;d<getDimensions();d++)
                    {
                        ASSERT_NEAR( data(i,d,t) , data2(i,d,t) ,1e-3);
                        ASSERT_GT( std::abs(data(i,d,t)) , 0 );
                    }
                }
        }

    }

    ASSERT_GT(success,0);
    EXPECT_LT(success,nSteps);

 }


TEST(configurations, IO)
{
    const int N = 1000;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();

    configurations.join(100,120);
    configurations.join(120,100);

    configurations.setHead(10,M+1);

    std::string filename {"testConf.h5"} ;

    configurations.saveHDF5(filename);



     auto configurations2 = pimc::pimcConfigurations::loadHDF5(filename); 
    
    auto & data2 = configurations2.dataTensor();


    ASSERT_EQ(configurations.nBeads() , configurations2.nBeads() );
    ASSERT_EQ(configurations.nParticles() , configurations2.nParticles() );

    for(int t=0;t<M+1;t++)
        for (int i=0;i<N;i++)
            for(int d=0;d<getDimensions();d++)
                {
                   ASSERT_NEAR( data(i,d,t) , data2(i,d,t), 1e-5);
                } 

                
    for(int i=0;i<N;i++)
    {
        ASSERT_EQ( configurations2.getChain(i).prev , configurations.getChain(i).prev ); 
        ASSERT_EQ( configurations2.getChain(i).next , configurations.getChain(i).next ); 
        ASSERT_EQ( configurations2.getChain(i).head , configurations.getChain(i).head ); 
        ASSERT_EQ( configurations2.getChain(i).tail , configurations.getChain(i).tail );
        ASSERT_EQ( configurations2.getGroups()[0].tails[0] , configurations2.getGroups()[0].tails[0] );
    } 

    configurations.saveHDF5(filename);

}

TEST(observables, IO)
{
    const int N = 1000;
    const int M = 50;   
    Real Beta=1;
    Real timeStep=Beta/M;

    pimc::geometryPBC_PIMC geo(300,300,300);



    pimc::particleGroup groupA{ 0 , N-1, N , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();

    data.setRandom();

    configurations.fillHeads();

    
     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    #if DIMENSIONS == 1
    auto V = pimc::makePotentialFunctor(
         [](Real x) {return 0.5*(x*x ) ;} ,
         [](Real x) {return 0.5*x ;} 
         );    
    #endif

    #if DIMENSIONS == 3

    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r ;} 
         );


    #endif


     std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);

    pimc::firstOrderAction S(sT,  sV);


    

    std::vector<std::shared_ptr<pimc::observable> > Os;
    {
        auto therm = std::make_shared<pimc::thermodynamicEnergyEstimator>();
        Os.push_back( std::make_shared<pimc::scalarObservable>(therm,std::string("eT") ) );

    }


    for (auto & O : Os)
    {
        O->accumulate(configurations,S);
        O->out(0);
        O->clear();
        O->out(1);
    }


    
}


TEST(configurations, worms)
{
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N + 1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});


    auto & data = configurations.dataTensor();
    const auto & tags = configurations.getTags();

    data.setRandom();

    configurations.fillHeads();

    configurations.save("randomConfiguration");


    int time=10;
    int iChain=0;

    // test open
    configurations.setHead( iChain , M );

    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );
    
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0],iChain);
    ASSERT_EQ(configurations.getGroupByChain(0).tails[0],iChain);


    ASSERT_TRUE( configurations.getChain(iChain).hasHead()     );
    ASSERT_TRUE(  configurations.getChain(iChain).hasTail() );

    // test close

    configurations.join(iChain,iChain);
    ASSERT_FALSE(configurations.getChain(iChain).hasHead() );
    ASSERT_FALSE(configurations.getChain(iChain).hasTail() );
    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 0 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 0 );

    // test creation of a new head and join
     iChain=5;
    configurations.setHead(iChain,M);
    configurations.pushChain(0);
    configurations.setHead(N,time);
    configurations.setTail(N,-1);
    TESTCHAIN(configurations, N, time, -1);
    TESTCHAIN(configurations, iChain, M, -1); 
    configurations.join(iChain,N);

    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );


    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , N );
    


    // test creation of a new tail and join

    iChain=5;
    configurations.setTail(iChain,-1);
    int newChain=configurations.pushChain(0);
    configurations.setHead(newChain,M);
    configurations.setTail(newChain,time);
    TESTCHAIN(configurations, newChain, M, time);
    TESTCHAIN(configurations, iChain, M, -1); 
    configurations.join(newChain,iChain); 


    ASSERT_EQ( configurations.getGroupByChain(0).heads.size() , 1 );
    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 1 );

    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , newChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , N );
    ASSERT_EQ(newChain,N+1);
    ASSERT_FALSE(configurations.getChain(N).hasTail() );
    ASSERT_FALSE(configurations.getChain(N+1).hasHead() );


    // test remove
    configurations.setTail(N,-1);

    ASSERT_EQ( configurations.getGroupByChain(0).tails.size() , 2 );

    configurations.setHead(N+1,0);
    
    ASSERT_EQ( configurations.getGroupByChain(0).tails[1] , N );
    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , N + 1 );
    
    configurations.removeChain(N); 
    configurations.removeChain(N);

    ASSERT_EQ( configurations.getGroupByChain(0).tails[0] , iChain );
    ASSERT_EQ( configurations.getGroupByChain(0).heads[0] , iChain );

    pimc::particleGroup groupB{ 0 , N-1, N + 2 , 1.0};

    pimc::pimcConfigurations configurations2(M , getDimensions() , {groupB});
    
    configurations2.setEnsamble(pimc::ensamble_t::grandCanonical);
    int newHead=M/2;
    int newTail=M/2 ;
    int iOpen=N/2;

    configurations2.dataTensor().setRandom();
    configurations2.fillHeads();
    int iNext=configurations2.getChain(iOpen).next;

    configurations2.setHead(iOpen,newHead);
    auto iChainTail=configurations2.pushChain(0);
    configurations2.setTail(iChainTail,newTail);
    configurations2.join(iChainTail,iNext);

    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );
    ASSERT_EQ( configurations2.getChain(iChainTail).next , iNext );
    ASSERT_EQ( configurations2.getChain(iNext).prev , iChainTail );
    ASSERT_EQ( configurations2.getGroupByChain(0).tails[0] , iChainTail );
    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );


    TESTCHAIN(configurations2,iOpen,newHead,-1);
    TESTCHAIN(configurations2,iChainTail,M,newTail);

    configurations2.setHead(iChainTail,M);
    configurations2.removeChain(iOpen);

    TESTCHAIN(configurations2,iOpen,M,newTail);


    ASSERT_EQ( configurations2.getGroupByChain(0).heads[0] , iOpen );
    ASSERT_EQ( configurations2.getGroupByChain(0).tails[0] , iOpen );


    configurations2.setTail(iOpen,-1);
    iChainTail=configurations2.pushChain(0);
    configurations2.setHead(iChainTail,M);
    configurations2.setTail(iChainTail,newTail);
    configurations2.join(iChainTail,iOpen);

    TESTCHAIN(configurations2,iOpen,M,-1);
    TESTCHAIN(configurations2,iChainTail,M,newTail);

    {
    pimc::particleGroup groupB{ 0 , 1 , 1 , 1.0};

    pimc::pimcConfigurations configurations3(M , getDimensions() , {groupB});

    configurations3.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations3.setChemicalPotential({0});
    
    

    
     auto & data3 = configurations3.dataTensor();
    data3.setRandom();


    pimc::pimcConfigurations configurations3Old(configurations3);
    const auto & data3Old = configurations3Old.dataTensor();
    
    const auto & tags3 = configurations3.getTags();
    const auto & tags3Old = configurations3Old.getTags();

    configurations3.setHeadTail(1,M,-1);
    configurations3.removeChain(1);


    for (int t=0;t<=M;t++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            ASSERT_NEAR( data3(0,d,t) , data3Old(0,d,t) ,1e-6);
            ASSERT_EQ( tags3(0,t) , tags3Old(0,t));
        }
    }

    int tHead=3;
    int newChain = configurations3.pushChain(0);
    configurations.setHeadTail(newChain,M,tHead);
    configurations3.copyData({tHead+1,M}  , 0, newChain  );

    for (int t=tHead+1;t<M;t++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            ASSERT_NEAR( data3(newChain,d,t) , data3Old(0,d,t) ,1e-6);
            ASSERT_EQ( tags3(newChain,t) , 1);
        }

      

    }

    for(int d=0;d<getDimensions();d++)
        {
        ASSERT_NEAR( data3(newChain,d,M) , data3Old(0,d,M) ,1e-6);
        ASSERT_EQ( tags3(newChain,M) , 0);

        }


    



    }


      {
        pimc::particleGroup groupB{ 0 , 2 , 3 , 1.0};

        pimc::pimcConfigurations configurations3(M , getDimensions() , {groupB});

        int iHead = 3;

        auto & data3 = configurations3.dataTensor();

        

        data3.setRandom();
        configurations3.fillHeads();
        configurations3.setHeadTail(0,iHead,-1);
        configurations3.setEnsamble(pimc::ensamble_t::grandCanonical);
        configurations3.setChemicalPotential({0});
        


        pimc::pimcConfigurations configurations3Old(configurations3);

        

        const auto & tags3Old = configurations3Old.getTags();
        const auto & tags3 = configurations3.getTags();
        const auto & data3Old = configurations3Old.dataTensor();

        
        configurations3.swap(0,1);


        for (int t=0;t<=M;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                ASSERT_FLOAT_EQ( data3(0,d,t) , data3Old(1,d,t) );
                
            }
            ASSERT_FLOAT_EQ( tags3(0,t) , tags3Old(1,t) );

        }

        ASSERT_EQ(configurations3.getGroups()[0].heads[0],1);
        ASSERT_EQ(configurations3.getGroups()[0].tails[0],1); 

      }

}


TEST(action,oneBodyGrandCanonical)
{
    int N=100;
    int M=50;

    pimc::particleGroup group{ 0 , N-1, N + 2 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {group});

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

    Real timeStep=1e-1;

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



TEST(run,free_harmonic_oscillator)
{   
    int N=1;
    int M=10;
    Real Beta = 1;


    pimc::geometryPBC_PIMC geo(300,300,300);

    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    

    configurations.dataTensor().setRandom();

    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();

    pimc::levyReconstructor reconstructor(M);

    pimc::levyMove freeMoves(5,0);

    Real delta=0.1;

    pimc::translateMove translMove(delta,(M+1)*N,0);
    

    Real C = 1e-1;
    int l = 3;

    
    pimc::openMove openMove(C,0,l);
    pimc::closeMove closeMove(C,0,l);

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    pimc::swapMove swapMove(l,N,0);

    pimc::tableMoves table;

    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    //table.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    //table.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");


    //table.push_back(& openMove,0.2,pimc::sector_t::diagonal,"open");
    
    table.push_back(& closeMove,0.2,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    table.push_back(& swapMove,0.8,pimc::sector_t::offDiagonal,"swap");
    

    randomGenerator_t randG(368);

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);
    

     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );


    Real R0=0.1;
    Real V0=1;



  
      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-(r*r));} ,
         [=](Real r) {return -2*r*V0*exp(-r*r)  ;}
          );




    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,N,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};

    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    pimc::firstOrderAction S(sT,  sV);
    
    int nTimes = 1000;
    int success = 0;
    int subSteps=1000;
    int correlationSteps=10;

   
    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(N, M);

    Real e=0;
    Real e2=0;

    std::ofstream f;
    std::ofstream fV;


    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    configurations.fillHeads();


    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat");
    fV.open("energyVirial.dat");

    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0;
        int nMeasurements=0;

        for (int k=0;k< subSteps;k++)
        {
            
            for (int j=0;j<correlationSteps;j++)
            {
                bool accepted=table.attemptMove(configurations, S, randG);

                if (accepted)
                {success+=1;}
            }
            
            if (!configurations.isOpen() )
            {
                Real tmp=energyEstimator(configurations,S);
                Real tmp1=viriralEnergy(configurations,S);
            
                nMeasurements++;
            
                eStep+=tmp;
                eVirialStep+=tmp1;
            }
            else
            {
               
            }
            
        }

        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;

        table >> std::cout;

        configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
    }

    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}


void accumulateBeadPosition(int i,std::array<Real,getDimensions()> & x, std::array<Real,getDimensions()> & x2, const pimc::configurations_t & configurations, const pimc::firstOrderAction & S)
    {
        const auto & data = configurations.dataTensor();
        for(int d=0;d<getDimensions();d++)
            {
                    x[d]+=data( 0 ,d, i) ;
                    x2[d]+=std::pow(data(0,d,i),2);

                }
    }


Real accumulateX2( const pimc::configurations_t & configurations, const pimc::firstOrderAction & S)
{
    Real x2=0;
    const auto & data = configurations.dataTensor();
    const auto & group = configurations.getGroups()[0];    

    for (int t=0;t<configurations.nBeads();t++)
        for (int i=group.iStart;i<=group.iEnd;i++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    x2+=data( i ,d, t)*data(i,d,t);
                }

            }
    

    return x2/configurations.nBeads();
            
}

Real accumulateX2SingleDirection(int iChain, const pimc::configurations_t & configurations, int direction, bool &  isCyclic, int t0 , int t1)
    {
        const auto & data = configurations.dataTensor();
        Real l2=0;

        int iCurrentChain=iChain;

        if (iCurrentChain < 0)
        {
            return 0;
        }
        do 
        {
            const auto & chain = configurations.getChain(iCurrentChain);

            for (int t=std::max(chain.tail + 1,t0);t<=std::min(chain.head , t1 );t++)
            {
                Real prefactor = ( (t == chain.tail+1) or (t == chain.head) ) ? 0.5 : 1;
                for(int d=0;d<getDimensions();d++)
                {
                        l2+=prefactor*data( iCurrentChain ,d, t)*data( iCurrentChain ,d, t);
                }
            }

            if (direction == 1)
            {
                iCurrentChain=chain.next;
            }
            else if (direction == -1)
            {
                iCurrentChain=chain.prev;
            }

        }
        while ( (iCurrentChain!= -1) and (iCurrentChain != iChain) );


        isCyclic= iCurrentChain == - 1 ? false : true;
        
        return l2;

    }

Real accumulateAverageLengthSquareSingleDirection(int iChain, const pimc::configurations_t & configurations, int direction, bool &  isCyclic, int t0 , int t1)
    {
        const auto & data = configurations.dataTensor();
        Real l2=0;

        int iCurrentChain=iChain;

        if (iCurrentChain < 0)
        {
            return 0;
        }
        do 
        {
            const auto & chain = configurations.getChain(iCurrentChain);

            for (int t=std::max(chain.tail + 1,t0);t<=std::min(chain.head - 1, t1 );t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                        l2+=std::pow(data( iCurrentChain ,d, t+1) - data(iCurrentChain,d,t),2);
                }
            }

            if (direction == 1)
            {
                iCurrentChain=chain.next;
            }
            else if (direction == -1)
            {
                iCurrentChain=chain.prev;
            }

        }
        while ( (iCurrentChain!= -1) and (iCurrentChain != iChain) );


        isCyclic= iCurrentChain == - 1 ? false : true;
        
        return l2;

    }


Real accumulateAverageLengthSquare(int iChain, const pimc::configurations_t & configurations, int t0, int t1)
{
    bool isCyclic;

    Real l2=accumulateAverageLengthSquareSingleDirection(iChain, configurations,+1,isCyclic,t0,t1);

    if (not isCyclic)
    {
        
        l2+=accumulateAverageLengthSquareSingleDirection(configurations.getChain(iChain).prev, configurations,-1,isCyclic,t0,t1);

        assert(isCyclic == false);

    }

    return l2;
}

int getWormLength(const pimc::configurations_t & configurations, int iGroup)
{
    if (not configurations.isOpen(iGroup) )
    {
        throw std::runtime_error("getWormLength requires an open sector");
    }

    int l=0;
    int iChainHead=configurations.getGroups()[iGroup].heads[0];

    int iChain=iChainHead;
    while ( iChain != -1  ) 
    {
        l+=  configurations.getChain(iChain).head - configurations.getChain(iChain).tail - 1;

        iChain=configurations.getChain(iChain).prev;
    }


    return l;
}



Real accumulateX2(int iChain, const pimc::configurations_t & configurations, int t0, int t1)
{
    bool isCyclic;

    Real l2=accumulateX2SingleDirection(iChain, configurations,+1,isCyclic,t0,t1);

    if (not isCyclic)
    {
        
        l2+=accumulateX2SingleDirection(configurations.getChain(iChain).prev, configurations,-1,isCyclic,t0,t1);

        assert(isCyclic == false);

    }

    return l2;
}




Real accumulateAverageLengthSquare(int iChain, const pimc::configurations_t & configurations)
{
    int t0=0;
    int t1=configurations.nBeads()-1;



    return accumulateAverageLengthSquare( iChain, configurations,0,t1);
}





auto meanBeadFixedLengths(int iChainBegin , int iChainEnd,  int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
{
    const auto & data = configurations.dataTensor();



    std::array<Real,getDimensions()> meanExpected;
    
    Real D = 0.5;
    Real mass = 1;
    int M= configurations.nBeads();

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        int l1 = (i - t0) > 0 ? i - t0 : i - t0 + M;
        int l2 = (t1 - i ) > 0 ? t1 - i  : t1 - i  + M;


        meanExpected[d]=(data(iChainBegin,d,t0)/ l1 + 
        data(iChainEnd,d,t1)/l2 )/(1./l1 + 1./l2 );

    }

    return meanExpected;
}

auto varianceBeadFixedLengths(int iChain , int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
{

    std::array<Real,getDimensions()> varianceExpected;
    const auto & data = configurations.dataTensor();

    Real D = 0.5;
    Real mass = 1;

    int M= configurations.nBeads();

    std::array<Real,3> difference;
    for (int d=0;d<getDimensions();d++)
    {
        difference[d]= data(iChain,d,t0) - data(iChain,d,t1);

        int l1 = (i - t0) > 0 ? i - t0 : i - t0 + M;
        int l2 = (t1 - i ) > 0 ? t1 - i  : t1 - i  + M;


       varianceExpected[d]=1./(1./l1 + 1./l2 )* 2 * D * timeStep / mass;




    }

    return varianceExpected;
}

TEST_F(configurationsTest,openCloseGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int t0=80;
    int l =40;

    SetUp(1,100,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    addOpenCloseFixedSegment(C , l , t0);

    int t1=(t0 + l)%M;
    int i = (t0 + 20)%M;

    Real nTrials = 1000000;
    Real nBurns = 100;

    std::array<Real,getDimensions() > x{  DLIST(0,0,0) };
    std::array<Real ,getDimensions() > x2{ DLIST(0,0,0) };

    auto & data = configurations.dataTensor();
    

    pimc::configurations_t configurationsInitial=configurations;

    const auto & dataInitial = configurations.dataTensor();


    accumulate(nBurns,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){accumulateBeadPosition(i,x,x2,confs,S) ;}  );


    Real openRatio = nOpen/nMoves;
    
    auto  meanExpected = meanBeadFixedLengths(0,0,t0,t1,i,timeStep,configurations);

    auto  varianceExpected = varianceBeadFixedLengths(0,t0,t1,i,timeStep,configurations);
    

  
    Real D = 0.5;
    Real mass = 1;


    std::array<Real,getDimensions()> difference;
    for (int d=0;d<getDimensions();d++)
    {
          difference[d]= data(0,d,t0) - data(0,d,t1);
        x[d]/=nClosed;
        x2[d]/=nClosed;

    }

    Real expectedOpenRatio= 1/( 1 + exp(pimc::freeParticleLogProbability(difference,S.getTimeStep()*l,mass))/C );

    std::cout << "Open ratio error : " << std::abs(openRatio - expectedOpenRatio)/expectedOpenRatio <<  std::endl;

    EXPECT_NEAR((openRatio - expectedOpenRatio)/expectedOpenRatio , 0 , 1e-2);


    std::cout << "Mean error: " << std::abs((x[0]  -  meanExpected[0])/meanExpected[0])<< std::endl;
    EXPECT_NEAR( (x[0]  -  meanExpected[0])/meanExpected[0] , 0 , 1e-2);

    std::cout << "var error : " << std::abs(x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] << std::endl;

    EXPECT_NEAR(( x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] , 0 , 3e-2);

}


TEST_F(configurationsTest,openChain)
{
    Real C=1;
    int nBeads=10;

    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);



    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);

    configurations.fillHeads();

    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.9,pimc::sector_t::offDiagonal);


    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 1000000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,10,pimc::sector_t::offDiagonal);

    l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);

    EXPECT_NEAR(l2 , 3* 2, 2*l2Error);

}


TEST_F(configurationsTest,openChain_manipulations)
{
     Real C=1;
    int nBeads=10;

    SetUp(3,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetSeed( time( NULL));

    SetRandom();

    int tHead=4;
    int l=60;
    int lShort=5;


    configurations.setHeadTail(2,tHead,-1);
    configurations.setHeadTail(0,M,8);
    configurations.join(0,1);
    configurations.join(1,2);

    auto & Spot=S.getPotentialAction();


    Real deltaS=Spot.evaluate(configurations,{tHead,std::min(tHead + l ,M) -1},1);
    

    Real sum=0;
    int i=1;

    const auto & data = configurations.dataTensor();


    for (int t=tHead  ; t<= M ; t++)
    {
        for (int d=0;d<getDimensions();d++)
        {
            Real prefactor = (t == tHead) or (t == M) ? 0.5 : 1;
            sum+=prefactor*0.5*data(i,d,t)*data(i,d,t);
        }
    }

    ASSERT_NEAR(deltaS,sum*timeStep,1e-6);

    int iPartner=2;

    deltaS=Spot.evaluate(configurations,{0,tHead + l - M  - 1},iPartner);

    sum=0;
    for (int t=0 ; t<=1 ; t++)
    {
        for (int d=0;d<getDimensions();d++)
        {
            Real prefactor = 0.5;
            sum+=prefactor*0.5*data(iPartner,d,t)*data(iPartner,d,t);
        }
    }

    ASSERT_NEAR(deltaS,sum*timeStep,1e-6);

}

TEST_F(configurationsTest,openChain_swap)
{
    Real C=1;
    int nBeads=10;
    SetUp(3,nBeads,1,{ DLIST(3000 , 3000 , 3000) });

    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetSeed( time(NULL));
    SetRandom();

    int l = 3;
    
    pimc::levyMove levy(l,0);
    pimc::swapMove swap( 5 , 3 , 0);
    //swap.setFixedLength();


    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    pimc::translateMove translate(0.3, 6*M , 0 );

    
    tab.push_back(&levy,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.2,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.2,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal, "swap");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal, "translate");


    int iHead=4;
    int iTail=8;


    configurations.setHeadTail(0,M,iTail );
    configurations.setHeadTail(1,M,- 1 );
    configurations.setHeadTail(2,iHead,- 1 );
    
    configurations.join(0,1);
    configurations.join(1,2);
    

/* 
    configurations.setHeadTail(0,M,iTail );
    configurations.setHeadTail(1,M,- 1 );
    configurations.setHeadTail(1,iHead,- 1 );
    configurations.join(0,1);
    configurations.join(1,2); */


    configurations.fillHeads();

    resetCounters();

    int nTrials=10000;
    int nBlocks=1000000;
    std::ofstream l2SwappedOut,l2NonSwappedOut,swapRatioOut,x2SwappedOut,x2NonSwappedOut;

    l2SwappedOut.open("l2Swapped.dat");
    l2NonSwappedOut.open("l2NonSwapped.dat");
    x2NonSwappedOut.open("x2NonSwapped.dat");
    x2SwappedOut.open("x2Swapped.dat");
    swapRatioOut.open("ratio.dat");

    for ( int iBurn=0;iBurn<10000;iBurn++)
    {
        tab.attemptMove(configurations,S,randG);
    }


    for (int i=0;i<nBlocks;i++)
    {
        Real l2Swapped=0;
        Real l2NonSwapped=0;
        Real confNonSwapped=0;
        Real confSwapped=0;
        Real x2Swapped=0;
        Real x2NonSwapped=0;

        for (int k=0;k<nTrials;k++)
            {
                for ( int iDec=0;iDec<10;iDec++)
                {
                    tab.attemptMove(configurations,S,randG);
                }


                Real tmp=accumulateAverageLengthSquare(0,configurations);
                Real tmpX2=accumulateX2(0,configurations,0,M);

                int nChainsLinked=configurations.buildPolimerList(0).size();

                if ( nChainsLinked == 2   )
                {
                    confNonSwapped+=1;
                    l2NonSwapped+=tmp;
                    x2NonSwapped += tmpX2;
                }
                else if ( nChainsLinked == 3   )
                {
                    confSwapped+=1;
                    l2Swapped+=tmp;
                    x2Swapped += tmpX2;

                }

            }
        
        l2SwappedOut << i << " " << l2Swapped/confSwapped << std::endl;
        l2NonSwappedOut << i << " " << l2NonSwapped/confNonSwapped << std::endl ;
        x2NonSwappedOut << i << " " << x2NonSwapped/confNonSwapped << std::endl ;
        x2SwappedOut << i << " " << x2Swapped/confSwapped << std::endl ;

        tab >> std::cout ;

        swapRatioOut << i << " " <<  confSwapped/(confNonSwapped + confSwapped) << std::endl;

    }

}

TEST_F(configurationsTest,openChain_advanceRecedeHead)
{
    Real C=1;
    int nBeads=10;
    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetSeed( time(NULL) );
    SetRandom();

    int t0=3;
    int l = 3;
    
    int iHead=8;
    int iTail=3;

    pimc::translateMove translate(0.3, 200*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    tab.push_back(&levy,1,pimc::sector_t::diagonal,"levy");

    tab.push_back(&levy,1,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&translate,1,pimc::sector_t::offDiagonal,"translate");
    


    configurations.setHeadTail(0,iHead,iTail);

    configurations.fillHeads();

    
    resetCounters();

    int nTrials=10000;
    int nBlocks=1000000;
    std::ofstream l2ShortOut,l2LongOut,swapRatioOut,x2ShortOut,x2LongOut;


    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    swapRatioOut.open("ratio.dat");
    x2ShortOut.open("x2Short.dat");
    x2LongOut.open("x2Long.dat");
    
    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);

    pimc::advanceHead advance(l,0);
    pimc::recedeHead recede(l,0);


    advance.setFixedLength();
    recede.setFixedLength();


    for (int i=0;i<nBlocks;i++)
    {
        Real l2Short=0;
        Real l2Long=0;
        Real x2Short=0;
        Real x2Long=0;
      
        Real nLong=0;
        Real nShort=0;


        for (int i=0;i<nTrials;i++)
        {
            for (int k=0;k<10;k++)
            {
                tab.attemptMove(configurations,S,randG);

            if (configurations.getChain(0).head == iHead)
            {
                advance.attemptMove(configurations,S,randG);
            }
            else
            {
                recede.attemptMove(configurations,S,randG);
            }
            
            }
            Real tmp=accumulateAverageLengthSquare(0,configurations);

            Real tmpX2=accumulateX2(0,configurations,0,M);


            if (configurations.getChain(0).head == iHead)
            {
                l2Short+=tmp;
                nShort+=1;
                x2Short+=tmpX2;
            }
            else
            {
                l2Long+=tmp;
                nLong+=1;
                x2Long+=tmpX2;
            }
        
        }

        swapRatioOut << i << " " << nShort/(nShort + nLong) << std::endl;
        l2ShortOut << i << " " << l2Short/nShort << std::endl;
        l2LongOut << i << " " << l2Long/nLong << std::endl ;
        x2LongOut << i << " " << x2Long/nLong << std::endl ;
        x2ShortOut << i << " " << x2Short/nShort << std::endl ;

        tab >> std::cout;

        const auto & sectorTab = tab.getOpenSectorTable(0);

        auto totTrial = std::accumulate( sectorTab.nTrials().begin(), sectorTab.nTrials().end(), 0);

        std::cout << "------- Fraction of proposals ------" << std::endl;

        for (int iMove=0;iMove<sectorTab.size();iMove++)
        {
            std::cout << sectorTab.names()[iMove] << "\t" << sectorTab.nTrials()[iMove]/totTrial << std::endl;
        }

        tab.resetCounters();


    }

}

TEST_F(configurationsTest,openChain_freeAdvanceRecede)
{
    Real C=1;
    int nBeads=10;
    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetSeed( time(NULL) );
    SetRandom();

    int t0=3;
    int l = 3;
    
    int iHead=8;
    int iTail=3;

    pimc::translateMove translate(0.3, 200*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    tab.push_back(&levy,1,pimc::sector_t::diagonal,"levy");

    tab.push_back(&levy,1,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&translate,1,pimc::sector_t::offDiagonal,"translate");

    configurations.setHeadTail(0,iHead,iTail);

    configurations.fillHeads();

    
    resetCounters();

    int nTrials=10000;
    int nBlocks=1000000;
    std::ofstream lOut;



    lOut.open("l.dat");
    
    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);


    pimc::advanceHead advance(l,0);
    pimc::recedeHead recede(l,0);


    //advance.setFixedLength();
    //recede.setFixedLength();

    tab.push_back(&advance,1,pimc::sector_t::offDiagonal,"advance");
    tab.push_back(&recede,1,pimc::sector_t::offDiagonal,"recede");



    for (int i=0;i<nBlocks;i++)
    {
        Real nMeasures=0;
        Real length=0;


        for (int i=0;i<nTrials;i++)
        {
            for (int k=0;k<10;k++)
            {
                tab.attemptMove(configurations,S,randG);    
            }
            
            int iChainHead=configurations.getGroups()[0].heads[0];
            int iChainTail=configurations.getGroups()[0].tails[0];

            

            int lCurrent = configurations.nParticles()*M;
            lCurrent+=configurations.getChain(iChainHead).head - (configurations.getChain(iChainHead).tail + 1);
            
            if (iChainTail != iChainHead)
            {
                lCurrent+=configurations.getChain(iChainTail).head - (configurations.getChain(iChainTail).tail + 1);    
            }
 
            nMeasures+=1;
            length+=lCurrent; 


        
        }

        lOut << i << " " << length/nMeasures << std::endl ;


        tab >> std::cout;

        
        tab.resetCounters();


    }

}

TEST_F(configurationsTest,levyReconstructor)
{
    int nBeads=10;
    int l = 5;
    int t0=0;
    int nTrials=100000;
    int nBlocks = 1000;

    SetUp(1,nBeads,1);
    SetUpNonInteractingHarmonicAction();
    SetGrandCanonicalEnsamble(0);

    auto tailBoundary = pimc::chainBoundary::free;
    auto headBoundary=pimc::chainBoundary::free;


    pimc::levyReconstructor reconstructor(2*l); 


    SetSeed( time(NULL) );

    SetRandom();

    auto & data = configurations.dataTensor();

        std::array x0 { 0, 0, 0 };
        std::array x1 { 1, 1, 1 };

        for(int d=0;d<getDimensions();d++)
        {
            data(0,d,t0)=x0[d];
            data(0,d,t0+l)=x1[d];
        }
        configurations.setHeadTail(0,t0+l,t0-1);

        reconstructor.setReconstructorBoundaries(tailBoundary,headBoundary );
        reconstructor.setGaussianParticleSampling();

        reconstructor.setSigma(1);


        // test reconstruction with fixed ends
        std::ofstream x2Out;

        x2Out.open("x2.out");

        x2Out << std::setprecision(5) ;



        
        std::array< int, 2 > timeRange { t0 , t0 + l};

        for(int iBlock=0;iBlock < nBlocks ; iBlock++)
        {
            Real x2=0;
            for (int t=0;t<nTrials;t++)
            {
                reconstructor.apply(configurations, timeRange , 0 ,  S , randG);
                x2+=accumulateX2(0,configurations,t0,t0+l);
            }

            x2Out << iBlock << " "  <<   x2/nTrials << std::endl;
        }

        x2Out.close();

        // 2.1 predicted value
    

}






TEST_F(configurationsTest,closedChains)
{
    Real C=1e-3;
    int nBeads=10;
    SetUp(12,nBeads,1);
    SetGrandCanonicalEnsamble(0.8 * getDimensions() /2. );
    SetUpNonInteractingHarmonicAction();


    //SetUpFreeParticleAction();
    SetSeed( 18 );
    SetRandom();
    

    int t0=3;
    int l = int( 0.7*nBeads);
    int lShort=int( 0.3*nBeads);

    int iHead=8;
    int iTail=3;

    pimc::translateMove translate(0.1, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lShort );
    pimc::closeMove close(C, 0, lShort );
    
    pimc::createWorm createWorm(C, 0, lShort , 1 );
    pimc::deleteWorm removeWorm(C, 0, lShort , 1);



    //open.setStartingBead(t0);
    //open.setStartingChain(0);

    //close.setStartingBead(t0);
    //close.setStartingChain(0);

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);

    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);


    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.2,pimc::sector_t::diagonal,"translate");
    tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");
    tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");



    tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&levy,0.2,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.05,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.05,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");

    tab.push_back(&advanceTail,0.05,pimc::sector_t::offDiagonal,"advanceTail");
    tab.push_back(&recedeTail,0.05,pimc::sector_t::offDiagonal,"recedeTail");
    tab.push_back(&advanceHead,0.1,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.1,pimc::sector_t::offDiagonal,"recedeHead");
    tab.push_back(&removeWorm,0.1,pimc::sector_t::offDiagonal,"removeWorm");



    configurations.fillHeads();

    //configurations.setHeadTail(0,10,-1);

    


    
    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    std::ofstream NOut,lOut,ratioOut;

    NOut.open("N.dat");
    lOut.open("l.dat");
    ratioOut.open("ratio.dat");
    

    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);


 


    //advance.setFixedLength();
    //recede.setFixedLength();

     NOut << std::setprecision(12);


    for (int i=0;i<nBlocks;i++)
    {
        Real N=0;
        Real nClosed=0;
        Real nOpen=0;
        Real length=0;


        for (int i=0;i<nTrials;i++)
        {
            for (int k=0;k<1;k++)
            {
                tab.attemptMove(configurations,S,randG);
            }
            
            if ( not configurations.isOpen() )
            {
                const auto & group= configurations.getGroups()[0];
                
                N+= (group.iEnd - group.iStart + 1);
                nClosed+=1;
            }
            else
            {
               
 
            nOpen+=1;
            
            }

        
        }


        NOut << i << " " << N << " " << nClosed << std::endl ;
        //lOut << i << " " << length/nOpen << std::endl ;

        ratioOut << i << " " << nOpen/(nOpen + nClosed) << std::endl ;

        
        tab >> std::cout;

        
        tab.resetCounters();        

    }

}

TEST_F(configurationsTest, openCloseSemiGrandCanonical)
{
    Real C=1e-3;
    int nBeads=10;
    SetUp( {1,1} ,nBeads,1);
    SetSemiGrandCanonicalEnsamble(  0   );
    //SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();


    SetSeed( 18 );
    SetRandom();

    int iStartA=configurations.getGroups()[0].iStart;
    int iStartB=configurations.getGroups()[1].iStart;
    

    int t0=8;
    int l = int( 0.7*nBeads);
    int lShort=int( 0.3*nBeads);

   
    pimc::translateMove translateA(0.1, 4000*M , 0 );
    pimc::translateMove translateB(0.1, 4000*M , 1 );

    pimc::levyMove levyA(l,0);
    pimc::levyMove levyB(l,1);
    

    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);

    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);


    pimc::openMove openA(C, 0, lShort );
    pimc::closeMove closeA(C, 0, lShort );


    openA.setStartingBead(t0);
    openA.setStartingChain(iStartA);
    openA.setLengthCut(lShort);

    closeA.setStartingBead(t0);
    closeA.setStartingChain(iStartA);
    closeA.setLengthCut(lShort);



    tab.push_back(&levyA,0.8,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translateA,0.1,pimc::sector_t::diagonal,"translateA");
     tab.push_back(&openA,0.1,pimc::sector_t::diagonal,"openA");

    
    tab.push_back(&levyB,0.8,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translateB,0.2,pimc::sector_t::diagonal,"translateB");
    
   



    tab.push_back(&levyA,0.7,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translateA,0.1,pimc::sector_t::offDiagonal,"translateA");
    tab.push_back(&closeA,0.1,pimc::sector_t::offDiagonal,"closeA");

    
    tab.push_back(&levyB,0.8,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.5,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translateB,0.1,pimc::sector_t::offDiagonal,"translateB");

    

    configurations.fillHeads();
    

   /*  configurations.setHeadTail(iStartA,t0,-1);
    configurations.setHeadTail(iStartA + 1,M,t0 + lShort-1 );
    configurations.join(iStartA+1,iStartA);

    configurations.setHeadTail(iStartB,t0 + lShort,t0-1);    */ 


    
    resetCounters();

    int nTrials=1000000;
    int nBlocks=100000;
    std::ofstream l2ClosedOut,l2OpenOut,ratioOut;

    l2ClosedOut.open("l2Closed.dat");
    l2OpenOut.open("l2Open.dat");
    ratioOut.open("ratio.dat");


    l2ClosedOut << std::setprecision(12);
    l2OpenOut << std::setprecision(12);


    int iChainSample=iStartB + 1;
    Real nClosed=0;
    Real nOpen=0;
    Real l2Closed=0;
    Real l2Open=0;
    Real ratioOpen=0;
    Real nRatioOpen=0;
    

    for (int iBlock=0;iBlock<nBlocks;iBlock++)
    {
        
        while ( (nClosed < nTrials) and (nOpen < nTrials) and ( nRatioOpen  < nTrials ) )
        {
            for (int k=0;k<1;k++)
            {
                tab.attemptMove(configurations,S,randG);
            }
            
            if (not ( ( configurations.getGroups()[1].iEnd < iChainSample ) and (iChainSample > 0) ) )
            {

                
                Real tmp=accumulateAverageLengthSquare(iChainSample,configurations);

                if ( not configurations.isOpen() )
                {
                    const auto & group= configurations.getGroups()[0];

                    
                    
                    l2Closed+=tmp;
                    nClosed+=1;
                }
                else
                {
                    l2Open+=tmp;
                    nOpen+=1;
                }
            }

            {
                if (  configurations.isOpen() )
                {
                    
                    ratioOpen+=1;
                }
                nRatioOpen+=1;
            }



        }


        if (nClosed == nTrials)
        {
            l2ClosedOut << iBlock << " " << l2Closed/nTrials << std::endl << std::flush;
            l2Closed=0;
            nClosed=0;
        }

        if (nOpen == nTrials)
        {
            l2OpenOut << iBlock << " " <<  l2Open/nTrials << std::endl << std::flush;
            l2Open=0;
            nOpen=0;
        }


        if ( nRatioOpen == nTrials )
        {
            ratioOut << iBlock << " " <<  ratioOpen/nTrials << std::endl << std::flush;
            ratioOpen=0;
            nRatioOpen=0;
        }


        tab >> std::cout;
        
        tab.resetCounters();

    }

}


TEST_F(configurationsTest, harmonicOscillatorMixture)
{
    Real C=1e-5;
    int nBeads=10;
    int N [] = {2 , 2};
    SetUp( {N[0],N[1]} ,nBeads,1);

    int buffer=1000;
    int nChains=(N[0] + buffer)*2;

    SetSemiGrandCanonicalEnsamble(  0   );
    //SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();


    SetSeed( time(NULL) );
    SetRandom();

    int iStartA=configurations.getGroups()[0].iStart;
    int iStartB=configurations.getGroups()[1].iStart;


    int t0=4;
    int l = int( 0.7*nBeads);
    int lShort=int( 0.3*nBeads);

    int iHead=8;
    int iTail=3;

    pimc::translateMove translateA(0.1, 4000*M , 0 );
    pimc::translateMove translateB(0.1, 4000*M , 1 );

    pimc::levyMove levyA(l,0);
    pimc::levyMove levyB(l,1);
    

    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);

    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);


    pimc::openMove openA(C, 0, lShort );
    pimc::closeMove closeA(C, 0, lShort );

    pimc::openMove openB(C, 1, lShort );
    pimc::closeMove closeB(C, 1, lShort );


    pimc::advanceHead advanceA(l,0);
    pimc::recedeHead recedeA(l,0);

    pimc::advanceHead advanceB(l,1);
    pimc::recedeHead recedeB(l,1);
    
    pimc::swapMove swapA( lShort , N[0] + 2 , 0);
    pimc::swapMove swapB( lShort , N[1] + 2, 1);


    //openA.setStartingBead(t0);
    //openA.setStartingChain(iStartA);
    //openA.setLengthCut(lShort);

    //closeA.setStartingBead(t0);
    //closeA.setStartingChain(iStartA);
    //closeA.setLengthCut(lShort);

    //closeB.setStartingBead(t0);
    //closeB.setStartingChain(iStartB);
    //closeB.setLengthCut(lShort);




    tab.push_back(&levyA,0.8,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translateA,0.1,pimc::sector_t::diagonal,"translateA");
    tab.push_back(&openA,0.1,pimc::sector_t::diagonal,"openA");

    tab.push_back(&levyB,0.9,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translateB,0.1,pimc::sector_t::diagonal,"translateB");
    //tab.push_back(&openB,0.1,pimc::sector_t::diagonal,"openB");

    tab.push_back(&levyA,0.6,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translateA,0.1,pimc::sector_t::offDiagonal,"translateA");
    tab.push_back(&swapA,0.1,pimc::sector_t::offDiagonal,"swapA");
    tab.push_back(&closeA,0.1,pimc::sector_t::offDiagonal,"closeA");


    
    tab.push_back(&levyB,0.7,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.05,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translateB,0.1,pimc::sector_t::offDiagonal,"translateB");
    //tab.push_back(&closeB,0.1,pimc::sector_t::offDiagonal,"closeB");
    tab.push_back(&swapB,0.1,pimc::sector_t::offDiagonal,"swapB");
    
    configurations.fillHeads();
    

   /*  configurations.setHeadTail(iStartA,t0,-1);
    configurations.setHeadTail(iStartA + 1,M,t0 + lShort-1 );
    configurations.join(iStartA+1,iStartA);

    configurations.setHeadTail(iStartB,t0 + lShort,t0-1);    */ 

    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(nChains, M);

    pimc::particleNumberEstimator particleEstimator0(0);
    pimc::particleNumberEstimator particleEstimator1(1);
    
    std::ofstream f;
    std::ofstream fV;
    std::ofstream fN;
    std::ofstream ratioOut;
    

    f.open("energy.dat");
    fV.open("energyVirial.dat");
    fN.open("N.dat");
    ratioOut.open("ratio.dat");

    

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    
    Real nClosed=0;
    Real nOpen=0;
    Real ratioOpen=0;
    Real nRatioOpen=0;

    Real e=0,eV=0,N0=0,N1=0;

    for (int iBlock=0;iBlock<nBlocks;iBlock++)
    {
        
        while ( (nClosed < nTrials) and (nOpen < nTrials) and ( nRatioOpen  < nTrials ) )
        {
            for (int k=0;k<1;k++)
            {
                tab.attemptMove(configurations,S,randG);
            }
            
            if ( not configurations.isOpen() )
            {

                e+=energyEstimator(configurations,S);
                eV+=viriralEnergy(configurations,S);

                N0+=particleEstimator0(configurations,S);
                N1+=particleEstimator1(configurations,S);
                

                nClosed+=1;
            }
            else
            {
                ratioOpen+=1;
                nOpen+=1;
            }


            nRatioOpen+=1;
        }

        if (nClosed == nTrials)
        {
        
            f << iBlock << "\t" << e/nTrials << std::endl << std::flush ;

            fV << iBlock << "\t" << eV/nTrials << std::endl << std::flush ;

            fN << iBlock << "\t" << N0/nTrials << "\t" << N1/nTrials << std::endl << std::flush;

            e=0;eV=0;N0=0;N1=0;
            nClosed=0;

        }

        if (nOpen == nTrials)
        {
            nOpen=0;

        }


        if ( nRatioOpen == nTrials )
        {
            ratioOut << iBlock << " " <<  ratioOpen/nTrials << std::endl << std::flush;
            ratioOpen=0;
            nRatioOpen=0;
        }


        tab >> std::cout;
        
        tab.resetCounters();

    }

}

TEST_F(configurationsTest, advanceRecedeSemiGrandCanonical)
{
    Real C=1e-3;
    int nBeads=10;
    SetUp( {2,2} ,nBeads,1);
    SetSemiGrandCanonicalEnsamble(  0   );
    //SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetSeed( 18 );
    SetRandom();

    int iStartA=configurations.getGroups()[0].iStart;
    int iStartB=configurations.getGroups()[1].iStart;

    int t0=9;
    int l = int( 0.7*nBeads);
    int lShort=int( 0.3*nBeads);

    int tHead=5 + 2;
    int tTail=3 + 2 ;
    

    pimc::translateMove translateA(0.1, 4000*M , 0 );
    pimc::translateMove translateB(0.1, 4000*M , 1 );

    pimc::levyMove levyA(l,0);
    pimc::levyMove levyB(l,1);
    

    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);

    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);




    tab.push_back(&levyA,0.7,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translateA,0.1,pimc::sector_t::offDiagonal,"translateA");

    tab.push_back(&levyB,0.8,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.5,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translateB,0.1,pimc::sector_t::offDiagonal,"translateB");

    configurations.fillHeads();
    
   

    configurations.setHeadTail( iStartA ,M, tTail);
    configurations.setHeadTail( iStartA + 1 ,tHead, -1);
    configurations.join(iStartA + 0,iStartA + 1);


    configurations.setHeadTail( iStartB , tTail + 1, -1);
    configurations.setHeadTail( iStartB + 1 , M  , tHead - 1);

    configurations.join(iStartB + 1,iStartB + 0);

    resetCounters();

    int nTrials=1000000;
    int nBlocks=100000;
    std::ofstream l2ShortOut,l2LongOut,ratioOut;

    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    ratioOut.open("ratio.dat");


    l2ShortOut << std::setprecision(12);
    l2LongOut << std::setprecision(12);

    int iChainSample=iStartB;
    Real nShort=0;
    Real nLong=0;
    Real l2Short=0;
    Real l2Long=0;
    Real ratioShort=0;
    Real nRatioShort=0;

    advanceHead.setFixedLength() ;
    recedeHead.setFixedLength() ;


    for (int iBlock=0;iBlock<nBlocks;iBlock++)
    {
        
        while ( (nShort < nTrials) and (nLong < nTrials) and ( nRatioShort  < nTrials ) )
        {
            for (int k=0;k<1;k++)
            {
                tab.attemptMove(configurations,S,randG);

                if (configurations.getChain(iStartA + 1).head  == tHead)
                {
                    advanceHead.attemptMove(configurations,S,randG);
                }
                else
                {
                    recedeHead.attemptMove(configurations,S,randG);
                }

            }
            
            
                
                Real tmp=accumulateAverageLengthSquare(iChainSample,configurations);

                if (configurations.getChain(iStartA + 1).head  == tHead)
                {
                    
                    
                    l2Short+=tmp;
                    nShort+=1;
                    ratioShort+=1;
                }
                else
                {
                    l2Long+=tmp;
                    nLong+=1;
                }
            

            

              nRatioShort+=1;


        }


        if (nShort == nTrials)
        {
            l2ShortOut << iBlock << " " << l2Short/nTrials << std::endl << std::flush;
            l2Short=0;
            nShort=0;
        }

        if (nLong == nTrials)
        {
            l2LongOut << iBlock << " " <<  l2Long/nTrials << std::endl << std::flush;
            l2Long=0;
            nLong=0;
        }


        if ( nRatioShort == nTrials )
        {
            ratioOut << iBlock << " " <<  ratioShort/nTrials << std::endl << std::flush;
            ratioShort=0;
            nRatioShort=0;
        }

        tab >> std::cout;

        tab.resetCounters();

    }
}

#include "../pimc/nConnectedChains.h"


TEST_F(configurationsTest,closedChain_twoBody)
{
    Real C=1e-3;
    int nBeads=10;
    int N=2;
    Real beta=0.1* nBeads;

    SetUp(N,nBeads,beta, { 300000} );

    //SetUpFreeParticleAction();    
    //SetUpNonInteractingHarmonicAction();
    
    //SetUpTwoBodyInteractionHarmonic();
    SetUpTwoBodyInteractionHarmonicInTrap();

    SetGrandCanonicalEnsamble(0.1 );
    SetSeed( time(NULL) );
    SetRandom();


   /*  auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return 0.*(r*r);} ,
         [=](Real r) {return 0*r  ;}
          );
    

    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    auto sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V2 ,geo,0,0);

    //S=pimc::firstOrderAction(sT,  sV2B); */

    int t0=7;
    int l = int( 0.8* 10);
    int lShort=int( 0.6* 10);
    int lOpen=lShort;

    pimc::translateMove translate(0.1, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lOpen );
    pimc::closeMove close(C, 0, lOpen );
    
    pimc::createWorm createWorm(C, 0, lShort , 1 );
    pimc::deleteWorm removeWorm(C, 0, lShort , 1);





    //open.setStartingBead(3);
    //open.setStartingChain(0);

    //close.setStartingBead(3);
    //close.setStartingChain(0);

    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);

    pimc::advanceHead advanceHead(1,0);
    pimc::recedeHead recedeHead(1,0);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();



    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);

    advanceHead.setMaximumParticleNumber(3);
    //recedeHead.setMinParticleNumber(1);

    pimc::nConnectedChains nConnectedChains;

    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.3,pimc::sector_t::diagonal,"translate");
    tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");
    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");
    
    tab.push_back(&levy,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&advanceHead,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.05,pimc::sector_t::offDiagonal,"recedeHead");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");


/*
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");
    tab.push_back(&advanceTail,0.05,pimc::sector_t::offDiagonal,"advanceTail");
    tab.push_back(&recedeTail,0.05,pimc::sector_t::offDiagonal,"recedeTail");
    tab.push_back(&advanceHead,0.1,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.1,pimc::sector_t::offDiagonal,"recedeHead");
    tab.push_back(&removeWorm,0.1,pimc::sector_t::offDiagonal,"removeWorm");
 */


    //int iHead = 7;
    //int lWormShort=10 + t0 - 3 ;

    //configurations.join(0,1);
    //configurations.join(1,0);

    //configurations.join(2,2);

    

    //configurations.join(1,0);    
    //configurations.setHead(1,0);

    configurations.setHeadTail(0,9,4);


    //configurations.setHeadTail(0,M,4-1);
    //configurations.join(0,1);
    //configurations.join(2,0);

    //configurations.join(1,1);

    //configurations.join(1,2);
    
    //configurations.join(0,1);
    //configurations.join(1,0);


    configurations.fillHeads();

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    
    std::ofstream NOut,l2ShortOut,l2LongOut, ratioOut,particleDistributionOut,wormDistributionOut,lWormOut,nConnectedChainsOut,nRingsOut;

    NOut.open("N.dat");
    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    particleDistributionOut.open("particleDistribution.dat");
    wormDistributionOut.open("wormDistribution.dat");
    lWormOut.open("lWorm.dat");
    nConnectedChainsOut.open("nConnected.dat");
    nRingsOut.open("nRings.dat");
    ratioOut.open("ratio.dat");

    

    std::vector<int > particleDistribution;
    int nMax=40;
    particleDistribution.resize(nMax,0);

    std::vector<int > wormDistribution;
    int nBeadsWormMax=400;
    wormDistribution.resize(nBeadsWormMax,0);


    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);

    NOut << std::setprecision(12);

    Real nShort=0;
    Real nLong=0;
    Real n=0;

    Real l2Long=0;
    Real l2Short=0;
    Real lWorm=0;

    Real nEstimator=0;
    Real r=0;

    Real nConnectedRingsEstimator=0;
    Real nConnectedChainsEstimator=0;

    for (int i=0;i<nBlocks;i++)
    {

        while (nShort<nTrials and nLong<nTrials and n<nTrials)
        {

            tab.attemptMove(configurations,S,randG);

           
           
            if ( configurations.isOpen() )
            //if (getWormLength(configurations,0) == lWormShort )
            {
               

                int l=getWormLength( configurations, 0);

                lWorm+=l;
                l2Short+=accumulateAverageLengthSquare(0,configurations);
                nShort+=1;
                r+=1;

                if (l < nBeadsWormMax)
                {
                    wormDistribution[l]+=1;
                }

                nConnectedChainsEstimator+=nConnectedChains.count(configurations,0);

            }
            else 
            {
                //l2Long+=accumulateAverageLengthSquare( 0,configurations );
                l2Long+=accumulateLengthSquare( configurations , {0,N-1} , {0,nBeads-1} , geo );
                nLong+=1;
                int currentN=configurations.nParticles();
                nEstimator+=currentN;
                if (currentN < nMax )
                {
                    particleDistribution[ currentN ] +=1 ;
                }

                nConnectedRingsEstimator+=nConnectedChains.count(configurations,0); 
            }

            n+=1;

        }

       if ( nShort == nTrials  )
       {
            l2ShortOut << i << " " << l2Short/nShort << std::endl << std::flush ;
            lWormOut << i << " " << lWorm/nShort << std::endl << std::flush;
            
            for (int iN=0;iN<nBeadsWormMax;iN++)
           {
               wormDistributionOut << i << " " << iN << " " <<  wormDistribution[iN] * 1./nTrials << " " << std::endl << std::flush ;
           }


           std::fill(wormDistribution.begin(), wormDistribution.end(), 0);

            nConnectedChainsOut << i << " " << nConnectedChainsEstimator/nTrials << std::endl << std::flush ;

            l2Short=0;
            nShort=0;
            lWorm=0;
            nConnectedChainsEstimator=0;

       }

       if (nLong == nTrials)
       {
           l2LongOut << i << " " << l2Long/nTrials << std::endl << std::flush;
           NOut << i << " " << nEstimator * 1./nTrials << std::endl << std::flush ;

           for (int iN=0;iN<nMax;iN++)
           {
               particleDistributionOut << i << " " << iN << " " <<  particleDistribution[iN] * 1./nTrials << " " << std::endl << std::flush ;
           }

           nRingsOut << i << " " << nConnectedRingsEstimator/nTrials << std::endl << std::flush ;



           std::fill(particleDistribution.begin(), particleDistribution.end(), 0);

           l2Long=0;
           nLong=0;
           nEstimator=0;
           nConnectedRingsEstimator=0;
            
       }

        if (n == nTrials)
        {
            ratioOut << i << " " << r/(n) << std::endl << std::flush ;
            r=0;
            n=0;
        }

        tab >> std::cout;
        tab.resetCounters();

    }

}

TEST_F(configurationsTest,createRemoveWorm)
{
    Real C=1;
    int nBeads=10;
    int N=1;
    Real beta=0.1* nBeads;
    
    SetUp(N,nBeads,beta, { 300000} );

    //SetUpFreeParticleAction();    
    SetUpNonInteractingHarmonicAction();
    
    //SetUpTwoBodyInteractionHarmonic();
    //SetUpTwoBodyInteractionHarmonicInTrap();

    SetGrandCanonicalEnsamble(0 );
    SetSeed( time(NULL) );
    SetRandom();

    
   /*  auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return 0.*(r*r);} ,
         [=](Real r) {return 0*r  ;}
          );
    

    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    auto sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V2 ,geo,0,0);

    //S=pimc::firstOrderAction(sT,  sV2B); */

    int t0=7;
    int l = int( 0.8* 10);
    int lShort=int( 0.6* 10);
    int lOpen=lShort;

    pimc::translateMove translate(0.1, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lOpen );
    pimc::closeMove close(C, 0, lOpen );
    
    pimc::createWorm createWorm(C, 0, 5 , 1 );
    pimc::deleteWorm removeWorm(C, 0, 5 , 1);

    createWorm.setStartingBead(5);
    createWorm.setFixedLength();
    createWorm.setGaussianSampling();

    removeWorm.setFixedLength();
    removeWorm.setStartingBead(5);
    removeWorm.setGaussianSampling();

    //open.setStartingBead(3);
    //open.setStartingChain(0);

    //close.setStartingBead(3);
    //close.setStartingChain(0);

    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);

    pimc::advanceHead advanceHead(1,0);
    pimc::recedeHead recedeHead(1,0);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();



    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);

    //advanceHead.setMaximumParticleNumber(3);
    //recedeHead.setMinParticleNumber(1);

    pimc::nConnectedChains nConnectedChains;

    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.3,pimc::sector_t::diagonal,"translate");
    //tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");
    tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");

    tab.push_back(&levy,0.6,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    //tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    //tab.push_back(&advanceHead,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    //tab.push_back(&recedeHead,0.05,pimc::sector_t::offDiagonal,"recedeHead");
    //tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");
    tab.push_back(&removeWorm,0.1,pimc::sector_t::offDiagonal,"removeWorm");

    //int iHead = 7;
    //int lWormShort=10 + t0 - 3 ;

    //configurations.join(0,1);
    //configurations.join(1,0);

    //configurations.join(2,2);

    

    //configurations.join(1,0);    
    //configurations.setHead(1,0);

    //configurations.setHeadTail(0,9,4);


    //configurations.setHeadTail(0,M,4-1);
    //configurations.join(0,1);
    //configurations.join(2,0);

    //configurations.join(1,1);

    //configurations.join(1,2);
    
    //configurations.join(0,1);
    //configurations.join(1,0);


    configurations.fillHeads();

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    
    std::ofstream NOut,l2ShortOut,l2LongOut, ratioOut,particleDistributionOut,wormDistributionOut,lWormOut,nConnectedChainsOut,nRingsOut;

    NOut.open("N.dat");
    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    particleDistributionOut.open("particleDistribution.dat");
    wormDistributionOut.open("wormDistribution.dat");
    lWormOut.open("lWorm.dat");
    nConnectedChainsOut.open("nConnected.dat");
    nRingsOut.open("nRings.dat");
    ratioOut.open("ratio.dat");

    

    std::vector<int > particleDistribution;
    int nMax=40;
    particleDistribution.resize(nMax,0);

    std::vector<int > wormDistribution;
    int nBeadsWormMax=400;
    wormDistribution.resize(nBeadsWormMax,0);


    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);

    NOut << std::setprecision(12);

    Real nShort=0;
    Real nLong=0;
    Real n=0;

    Real l2Long=0;
    Real l2Short=0;
    Real lWorm=0;

    Real nEstimator=0;
    Real r=0;

    Real nConnectedRingsEstimator=0;
    Real nConnectedChainsEstimator=0;

    for (int i=0;i<nBlocks;i++)
    {

        while (nShort<nTrials and nLong<nTrials and n<nTrials)
        {

            tab.attemptMove(configurations,S,randG);

            if ( configurations.isOpen() )
            //if (getWormLength(configurations,0) == lWormShort )
            {
               
                int iChainTail=configurations.getGroups()[0].tails[0];

                int l=getWormLength( configurations,0 );

                lWorm+=l;
                l2Short+=accumulateAverageLengthSquare(iChainTail,configurations);
                nShort+=1;
                r+=1;

                if (l < nBeadsWormMax)
                {
                    wormDistribution[l]+=1;
                }

                nConnectedChainsEstimator+=nConnectedChains.count(configurations,0);

            }
            else 
            {
                //l2Long+=accumulateAverageLengthSquare( 0,configurations );
                l2Long+=accumulateLengthSquare( configurations , {0,N-1} , {0,nBeads-1} , geo );
                nLong+=1;
                int currentN=configurations.nParticles();
                nEstimator+=currentN;
                if (currentN < nMax )
                {
                    particleDistribution[ currentN ] +=1 ;
                }

                nConnectedRingsEstimator+=nConnectedChains.count(configurations,0); 
            }

            n+=1;

        }

       if ( nShort == nTrials  )
       {
            l2ShortOut << i << " " << l2Short/nShort << std::endl << std::flush ;
            lWormOut << i << " " << lWorm/nShort << std::endl << std::flush;
            
            for (int iN=0;iN<nBeadsWormMax;iN++)
           {
               wormDistributionOut << i << " " << iN << " " <<  wormDistribution[iN] * 1./nTrials << " " << std::endl << std::flush ;
           }


           std::fill(wormDistribution.begin(), wormDistribution.end(), 0);

            nConnectedChainsOut << i << " " << nConnectedChainsEstimator/nTrials << std::endl << std::flush ;

            l2Short=0;
            nShort=0;
            lWorm=0;
            nConnectedChainsEstimator=0;

       }

       if (nLong == nTrials)
       {
           l2LongOut << i << " " << l2Long/nTrials << std::endl << std::flush;
           NOut << i << " " << nEstimator * 1./nTrials << std::endl << std::flush ;

           for (int iN=0;iN<nMax;iN++)
           {
               particleDistributionOut << i << " " << iN << " " <<  particleDistribution[iN] * 1./nTrials << " " << std::endl << std::flush ;
           }

           nRingsOut << i << " " << nConnectedRingsEstimator/nTrials << std::endl << std::flush ;



           std::fill(particleDistribution.begin(), particleDistribution.end(), 0);

           l2Long=0;
           nLong=0;
           nEstimator=0;
           nConnectedRingsEstimator=0;
            
       }

        if (n == nTrials)
        {
            ratioOut << i << " " << r/(n) << std::endl << std::flush ;
            r=0;
            n=0;
        }
        
        tab >> std::cout;

        
        tab.resetCounters();

    }

}

TEST_F(configurationsTest,swap_twoBody)
{
    Real C=1e-3;
    int nBeads=10;
    int N=2;
    Real beta=0.1* nBeads;

    SetUp(N,nBeads,beta, { 300000} );

    //SetUpFreeParticleAction();    
    SetUpNonInteractingHarmonicAction();

    //SetUpTwoBodyInteractionHarmonic();
    //SetUpTwoBodyInteractionHarmonicInTrap();
    
    SetGrandCanonicalEnsamble(0 );

    SetSeed( time(NULL) );
    SetRandom();



   /*  auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return 0.*(r*r);} ,
         [=](Real r) {return 0*r  ;}
          );
    

    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    auto sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V2 ,geo,0,0);

    //S=pimc::firstOrderAction(sT,  sV2B); */

    int l = int( 0.8* 10);
    int lShort=int( 0.4* 10);
    int lOpen=3;

    pimc::translateMove translate(0.1, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lOpen );
    pimc::closeMove close(C, 0, lOpen );
    
    pimc::createWorm createWorm(C, 0, lShort , 1 );
    pimc::deleteWorm removeWorm(C, 0, lShort , 1);

    //open.setStartingBead(t0);
    //open.setStartingChain(0);

    //close.setStartingBead(t0);
    //close.setStartingChain(0);

    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);

    pimc::advanceHead advanceHead(lOpen,0);
    pimc::recedeHead recedeHead(lOpen,0);

    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();



    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);
    swap.setFixedLength();



    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.3,pimc::sector_t::diagonal,"translate");
    tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");
    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");

    tab.push_back(&levy,0.6,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    //tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");


    int iHead = 10;
    //int lWormShort=10 + t0 - 3 ;

    //configurations.join(1,0);    
    //configurations.setHead(1,0);

    configurations.setHeadTail(0,iHead,-1);

    //configurations.join(0,1);

    //configurations.join(1,2);
    
    configurations.fillHeads();
    

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    std::ofstream l2ShortOut,l2LongOut,ratioOut;

    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    
    

    ratioOut.open("ratio.dat");

    std::vector<int > particleDistribution;
    int nMax=40;
    particleDistribution.resize(nMax,0);

    std::vector<int > wormDistribution;
    int nBeadsWormMax=400;
    wormDistribution.resize(nBeadsWormMax,0);


    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);


    Real nShort=0;
    Real nLong=0;
    Real n=0;

    Real l2Long=0;
    Real l2Short=0;

    Real r=0;

    for (int i=0;i<nBlocks;i++)
    {

        while (nShort<nTrials and nLong<nTrials and n<nTrials)
        {
            
            tab.attemptMove(configurations,S,randG);
    
            if ( configurations.getGroups()[0].heads[0] == 0 )
            {
                
                l2Short+=accumulateAverageLengthSquare(0,configurations);
                nShort+=1;
                r+=1;

                
                
            }
            else if ( configurations.getGroups()[0].heads[1] == 0 )
            {
                
                l2Long+=accumulateAverageLengthSquare( 0,configurations );
                nLong+=1;     

            } 

            n+=1;

        }

       if ( nShort == nTrials  )
       {
            l2ShortOut << i << " " << l2Short/nShort << std::endl << std::flush ;
            


           
            l2Short=0;
            nShort=0;
           
       }

       if (nLong == nTrials)
       {
           l2LongOut << i << " " << l2Long/nTrials << std::endl << std::flush;


           
           l2Long=0;
           nLong=0;
                  }

        if (n == nTrials)
        {
            ratioOut << i << " " << r/(n) << std::endl << std::flush ;
            r=0;
            n=0;
        }
        
        tab >> std::cout;

        
        tab.resetCounters();

    }

}

TEST_F(configurationsTest,driver)
{
    Real C=1e-3;
    int nBeads=10;
    SetUp(12,nBeads,1);
    SetGrandCanonicalEnsamble(0.8 * getDimensions() /2. );
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetSeed( 345 );
    SetRandom();

    std::ifstream ifs("/home/luca/source/qmc3/run/input.json");


    pimc::json_t j = pimc::json_t::parse(ifs);


    pimc::ensamble_t currentEnsamble = pimc::grandCanonical;
    std::vector<int> nMaxParticles;   
    if (currentEnsamble == pimc::grandCanonical)
    {
         
        nMaxParticles = j["maxParticles"].get<std::vector<int> >();

    }
    else
    {
        nMaxParticles = j["particles"].get<std::vector<int> >();
    }
    




    // ------------------ build moves ------------------------------



    pimc::moveConstructor pimcMoveConstructor(nMaxParticles,nBeads);
    pimcMoveConstructor.registerMove<pimc::levyMove>("levy");
    pimcMoveConstructor.registerMove<pimc::translateMove>("translate");
    pimcMoveConstructor.registerMove<pimc::openMove>("open");
    pimcMoveConstructor.registerMove<pimc::closeMove>("close");
    pimcMoveConstructor.registerMove<pimc::moveHead>("moveHead");
    pimcMoveConstructor.registerMove<pimc::moveTail>("moveTail");
    pimcMoveConstructor.registerMove<pimc::swapMove>("swap");

    if (currentEnsamble == pimc::ensamble_t::grandCanonical )
    {
        pimcMoveConstructor.registerMove<pimc::advanceHead>("advanceHead");
        pimcMoveConstructor.registerMove<pimc::recedeHead>("recedeHead");
        pimcMoveConstructor.registerMove<pimc::createWorm>("createWorm");
        pimcMoveConstructor.registerMove<pimc::deleteWorm>("deleteWorm");
    }




    int t0=3;
    int l = int( 0.7*nBeads);
    int lShort=int( 0.3*nBeads);


    int iHead=8;
    int iTail=3;

    pimc::translateMove translate(0.3, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lShort );
    pimc::closeMove close(C, 0, lShort );
    
    pimc::createWorm createWorm(C, 0, lShort , 1 );
    pimc::deleteWorm removeWorm(C, 0, lShort , 1);



    //open.setStartingBead(t0);
    //open.setStartingChain(0);

    //close.setStartingBead(t0);
    //close.setStartingChain(0);



     auto jAdvance = R"({
                "kind": "advanceHead",
                "reconstructionMaxLength": 3,
                "set" : 0
            })"_json;

    auto jRecede = R"({
                "kind": "recedeHead",
                "reconstructionMaxLength": 3,
                "set" : 0
            })"_json;
    

    auto jOpen= R"( 
         {
                "kind": "open",
                "reconstructionMaxLength": 3,
                "C": 0.001,
                "set" : 0,
                "particles" : [2000]
        }

    )"_json;

    auto jClose= R"( 
         {
                "kind": "close",
                "reconstructionMaxLength": 3,
                "C": 0.001,
                "set" : 0,
                "particles" : [2000]
        }

    )"_json;

    auto jCreateWorm= R"( 
        {
            "kind": "deleteWorm",
            "reconstructionMaxLength": 3,
            "C": 0.001,
            "alpha": 1,
            "set" : 0,
            "particles" : [2000]
            
        }

    )"_json;

    auto jDeleteWorm= R"( 
        {
            "kind": "deleteWorm",
            "reconstructionMaxLength": 3,
            "C": 0.001,
            "alpha": 1,
            "set" : 0,
            "particles" : [2000]
            
        }

    )"_json;

    auto jTranslate= R"( 
       {
        "kind": "translate",
        "delta": 0.1,
        "set" : 0,
        "particles" : [2000]    
        }

    )"_json;


    auto jSwap = R"({
                "kind": "swap",
                "reconstructionMaxLength": 3,
                "set" : 0,
                "particles" : [2000]
            })"_json;


    pimc::advanceHead advanceHead(jAdvance);

    pimc::recedeHead recedeHead(jRecede);

    
    pimc::swapMove swap( jSwap);




    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::diagonal,"translate");
    tab.push_back(&open,0.2,pimc::sector_t::diagonal,"open");
    tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");



    tab.push_back(&close,0.2,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&levy,0.1,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");

    //tab.push_back(&advanceTail,1,pimc::sector_t::offDiagonal,"advanceTail");
    //tab.push_back(&recedeTail,1,pimc::sector_t::offDiagonal,"recedeTail");
    tab.push_back(&advanceHead,0.1,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.1,pimc::sector_t::offDiagonal,"recedeHead");
    tab.push_back(&removeWorm,0.1,pimc::sector_t::offDiagonal,"removeWorm");


    configurations.fillHeads();

    //tab = pimcMoveConstructor.createTable( j["movesTable"] );

    
    resetCounters();

    int nTrials=100000;
    int nBlocks=1000000;
    std::ofstream NOut,lOut,ratioOut;

    NOut.open("N.dat");
    lOut.open("l.dat");
    ratioOut.open("ratio.dat");

    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);


    //advance.setFixedLength();
    //recede.setFixedLength();

    for (int i=0;i<nBlocks;i++)
    {
        Real N=0;
        Real nClosed=0;
        Real nOpen=0;
        Real length=0;

        for (int i=0;i<nTrials;i++)
        {
            for (int k=0;k<10;k++)
            {
                tab.attemptMove(configurations,S,randG);    
            }
            
            if ( not configurations.isOpen() )
            {
                const auto & group= configurations.getGroups()[0];

                N+= group.iEnd - group.iStart + 1;
                nClosed+=1;
            }
            else
            {
                int iChainHead=configurations.getGroups()[0].heads[0];
                int iChainTail=configurations.getGroups()[0].tails[0];


                int lCurrent = configurations.nParticles()*M;

                lCurrent+=configurations.getChain(iChainHead).head - (configurations.getChain(iChainHead).tail + 1);
            
                if (iChainTail != iChainHead)
                {
                    lCurrent+=configurations.getChain(iChainTail).head - (configurations.getChain(iChainTail).tail + 1);    
                }
 
            nOpen+=1;
            length+=lCurrent; 
            }
            
        
        }

        NOut << i << " " << N/nClosed << std::endl ;
        lOut << i << " " << length/nOpen << std::endl ;

        ratioOut << i << " " << nOpen/(nOpen + nClosed) << std::endl ;

        tab >> std::cout;
        
        tab.resetCounters();


    }

}



/* TEST_F(configurationsTest,openChain_advanceRecede_tail)
{
    Real C=1;
    int nBeads=10;
    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetSeed( time(NULL) );
    SetRandom();
    
    

    int t0=3;
    int l = 3;
    
    int iHead=5;
    int iTail=8;


    pimc::translateMove translate(0.3, 200*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    tab.push_back(&levy,1,pimc::sector_t::diagonal,"levy");

    tab.push_back(&levy,1,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&translate,1,pimc::sector_t::offDiagonal,"translate");
    


    configurations.setHeadTail(0,M,iTail);
    configurations.setHeadTail(1,iHead,-1);
    configurations.join(0,1);


    configurations.fillHeads();

    
    resetCounters();

    int nTrials=10000;
    int nBlocks=1000000;
    std::ofstream l2ShortOut,l2LongOut,swapRatioOut,x2ShortOut,x2LongOut;


    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    swapRatioOut.open("ratio.dat");
    x2ShortOut.open("x2Short.dat");
    x2LongOut.open("x2Long.dat");
    
    //pimc::advanceHeadTest advance(l);
    //pimc::recedeHeadTest recede(l);

    pimc::advanceTail advance(l,0);
    pimc::recedeTail recede(l,0);


    advance.setFixedLength();
    recede.setFixedLength();


    for (int i=0;i<nBlocks;i++)
    {
        Real l2Short=0;
        Real l2Long=0;
        Real x2Short=0;
        Real x2Long=0;
      
        Real nLong=0;
        Real nShort=0;


        for (int i=0;i<nTrials;i++)
        {
            for (int k=0;k<10;k++)
            {
                tab.attemptMove(configurations,S,randG);


            int iChainTail=configurations.getGroups()[0].tails[0];
            int iTailCurrent=configurations.getChain(iChainTail).tail;
            
            if ( iTailCurrent == iTail)
            {
                advance.attemptMove(configurations,S,randG);
            }
            else
            {
                recede.attemptMove(configurations,S,randG);
            }
            
            }

            int iChainTail=configurations.getGroups()[0].tails[0];
            int iTailCurrent=configurations.getChain(iChainTail).tail;


            Real tmp=accumulateAverageLengthSquare(iChainTail,configurations);

            Real tmpX2=accumulateX2(iChainTail,configurations,0,M);

           
           
            if (iTailCurrent == (iTail + l)%M )
            {
                l2Short+=tmp;
                nShort+=1;
                x2Short+=tmpX2;
            }
            else
            {
                l2Long+=tmp;
                nLong+=1;
                x2Long+=tmpX2;
            }
        
        }

        swapRatioOut << i << " " << nShort/(nShort + nLong) << std::endl;
        l2ShortOut << i << " " << l2Short/nShort << std::endl;
        l2LongOut << i << " " << l2Long/nLong << std::endl ;
        x2LongOut << i << " " << x2Long/nLong << std::endl ;
        x2ShortOut << i << " " << x2Short/nShort << std::endl ;

        tab >> std::cout;

        const auto & sectorTab = tab.getOpenSectorTable(0);

        auto totTrial = std::accumulate( sectorTab.nTrials().begin(), sectorTab.nTrials().end(), 0);

        //std::cout << configurations.getChain(0).tail << std::endl;

        tab.resetCounters();


    }

} */


TEST_F(configurationsTest,closedChain_free)
{
    Real C=1;
    int nBeads=10;

    SetUp(2,nBeads,1);
    //SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);




    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);
    configurations.join(1,0);


    //configurations.fillHeads();

    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);

    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 100000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::diagonal);

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    
    EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);


}


TEST_F(configurationsTest,closedChain_harmonic)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    //SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();

    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);

/* 

    configurations.setTail(0,-1);
    configurations.setHead(1,M);
    configurations.join(0,1);
    configurations.join(1,0);


    //configurations.fillHeads();
 */
    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);

    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;
    int i=4;

    int nTrials = 100000000;
    std::array<Real,getDimensions()> x2 { DLIST(0,0,0) };
    std::array<Real,getDimensions() > x { DLIST(0,0,0) };


    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;
    
    accumulateBeadPosition(i, x, x2, confs, S);
    
    }  ,100,pimc::sector_t::diagonal);

    for (int d=0;d<getDimensions();d++)
    {
        x[d]/=nClosed;
        x2[d]/=nClosed;
    }


    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    }

}

TEST_F(configurationsTest,openClosedChain_harmonic)
{
    Real C=1;
    int nBeads=10;
    int order =2;

    SetUp(1,nBeads,1, { DLIST(3000 , 3000 , 3000) });
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction(order);
    //SetUpFreeParticleAction();
    SetSeed(time(NULL));
    SetRandom();

    int t0=5;
    int l = 8;

    int lOpen=3;

    pimc::levyMove levy(l,0);

    assert(configurations.getChemicalPotential(0) == 0 );
    pimc::openMove open(C, 0, lOpen );
    pimc::closeMove close(C, 0, lOpen );

    //pimc::openMoveTest open(C, 0, lOpen,9 );
    //pimc::closeMoveTest close(C, 0, lOpen,9 );
    //close.setOrder(2);
    //open.setOrder(2);
    
    
    pimc::translateMove translate(0.3, 4*M , 0 );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    //open.setStartingBead(t0);
    //close.setStartingBead(t0);


    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);


    //configurations.join(0,2);
    //configurations.join(2,1);
    //configurations.join(1,3);
    //configurations.join(3,0);

    
    
    


    tab.push_back(&levy,0.7,pimc::sector_t::diagonal , "levy");
    tab.push_back(&translate,0.1,pimc::sector_t::diagonal,"translate");
    tab.push_back(&open,0.2,pimc::sector_t::diagonal,"open");
    
    tab.push_back(&translate,0.2,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&levy,0.2,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.2,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.2,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&close,0.2,pimc::sector_t::offDiagonal,"close");
    



    //configurations.setHeadTail(0,M,-1);

    //configurations.setTail(0,-1);
    //configurations.setHead(1,M);
    //configurations.join(0,1);
    //configurations.join(1,0);


    configurations.fillHeads();


    //Real l2=0;
    //Real l2Var=0;
    //Real l2Error=0;
   // int iNewHead=4;

    int nTrials = 1000000;

    std::ofstream l2OpenOut,l2ClosedOut,openFractionOut,x2OpenOut,x2ClosedOut;

    l2OpenOut.open("l2Open.dat");
    l2ClosedOut.open("l2Closed.dat");
    
    openFractionOut.open("openFraction.dat");
    x2OpenOut.open("x2Open.dat");
    x2ClosedOut.open("x2Closed.dat");

    pimc::sampleSingleClosedHarmonicChain sampleClosedChain( 1 ,0, M);

    

    pimc::sampleSingleOpenHarmonicChain sampleOpenChain( 1 ,0, M);

    for (int iBlock=0;iBlock < 10000000 ; iBlock++)
    {
        Real l2Closed=0;
        Real l2Open=0;
        Real x2Closed=0;
        Real x2Open=0;

        Real nOpen=0;
        Real nClosed=0;

        Real nClosedChainAttempts=0;
        Real nClosedChainSuccess=0;
        
        Real nOpenChainAttempts=0;
        Real nOpenChainSuccess=0;

        Real nCloseAttempts=0;
        Real nOpenAttempts=0;

        Real nOpenSuccess=0;
        Real nCloseSuccess=0;

        for (int n=0;n<nTrials;n++)
        {
            for (int k=0;k<1;k++)
            {
    

                tab.attemptMove(configurations,S, randG);

            }
            

            auto tmp=accumulateAverageLengthSquare(0,configurations) ;

            auto tmpX2=accumulateX2(0,configurations,0,M) ;


            if (! configurations.isOpen(0) )
            {
                l2Closed+=tmp;
                nClosed+=1;
                x2Closed+=tmpX2;
               
            }
            else
            {
                x2Open+=tmpX2;
                l2Open+=tmp;
                nOpen+=1;
            }
        
        } 

        l2ClosedOut << iBlock << "\t" << l2Closed/nClosed << std::endl;
        x2ClosedOut << iBlock << "\t" << x2Closed/nClosed << std::endl;
        x2OpenOut << iBlock << "\t" << x2Open/nOpen << std::endl;
        
        l2OpenOut << iBlock << "\t" << l2Open/nOpen << std::endl;
        openFractionOut << iBlock <<  "\t" << nOpen/(nOpen + nClosed) << std::endl;

        tab >> std::cout;

        if (nClosedChainAttempts > 0)
        {
            std::cout << "Acceptance closed sampling: " << nClosedChainSuccess/nClosedChainAttempts << std::endl;
        }
         if (nOpenChainAttempts > 0)
        {
            std::cout << "Acceptance open sampling: " << nOpenChainSuccess/nOpenChainAttempts << std::endl;
        }
        if (nOpenAttempts > 0)
        {
            std::cout << "Acceptance open->close:  " << nOpenSuccess/nOpenAttempts << std::endl;
        }
        if (nCloseAttempts > 0)
        {
            std::cout << "Acceptance closed-> open:  " << nCloseSuccess/nCloseAttempts << std::endl;
        }
    }

    x2OpenOut.close();
    x2ClosedOut.close();
    
    l2OpenOut.close();
    l2ClosedOut.close();
    openFractionOut.close();

}


TEST_F(configurationsTest,createRemove_harmonic)
{
    Real C=1;
    int nBeads=10;
    int order =2;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction(order);

    SetSeed(time(NULL));

    SetRandom();
    
    int t0=7;
    int l = 4;


    pimc::levyMove levy(l,0);

    pimc::createWorm open(C, 0, l , 1 );
    pimc::deleteWorm close(C, 0, l , 1);


    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setFixedLength();
    close.setFixedLength();


    pimc::translateMove translate(0.1, 3*M , 0 );


    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::diagonal,"translate");
     tab.push_back(&open,0.3,pimc::sector_t::diagonal,"create");
   


    tab.push_back(&levy,0.2,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.2,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.2,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&close,0.3,pimc::sector_t::offDiagonal,"remove");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");


    //configurations.setHeadTail(0,M,-1);

    //configurations.setTail(0,-1);
    //configurations.setHead(1,M);
    //configurations.join(0,1);
    //configurations.join(1,0);


    //configurations.fillHeads();


    //Real l2=0;
    //Real l2Var=0;
    //Real l2Error=0;
   // int iNewHead=4;

    int nTrials = 10000;

    std::ofstream l2OpenOut,l2ClosedOut,openFractionOut,x2OpenOut,x2ClosedOut;


    l2OpenOut.open("l2Open.dat");
    l2ClosedOut.open("l2Closed.dat");
    x2OpenOut.open("x2Open.dat");
    x2ClosedOut.open("x2Closed.dat");
    
    openFractionOut.open("openFraction.dat");


    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
        Real l2Closed=0;
        Real l2Open=0;
        nOpen=0;
        nClosed=0;
        Real x2Closed=0;
        Real x2Open=0;

        for (int n=0;n<nTrials;n++)
        {
            for (int k=0;k<50;k++)
            {
                bool accept = tab.attemptMove(configurations,S,randG);
            }

            if (! configurations.isOpen(0) )
            {
                auto tmp=accumulateAverageLengthSquare(0,configurations) ;
                x2Closed+=accumulateX2(0, configurations, 0, M);
                l2Closed+=tmp;
                nClosed+=1;
               
            }
            else
            {
                auto tmp=accumulateAverageLengthSquare(1,configurations) ;
                l2Open+=tmp;
                nOpen+=1;
                x2Open+=accumulateX2(1, configurations, 0, M);
            }
        
        } 

        l2ClosedOut << iBlock << "\t" << l2Closed/nClosed << std::endl;
        x2ClosedOut << iBlock << "\t" << x2Closed/nClosed << std::endl;
        x2OpenOut << iBlock << "\t" << x2Open/nOpen << std::endl;
        
        l2OpenOut << iBlock << "\t" << l2Open/nOpen << std::endl;
        openFractionOut << iBlock <<  "\t" << nOpen/(nOpen + nClosed) << std::endl;

        tab >> std::cout;

        
    }


    l2OpenOut.close();
    l2ClosedOut.close();
    openFractionOut.close();



}


TEST_F(configurationsTest,openClosedChain_free)
{
    Real C=1e-1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();

    SetSeed(time(NULL));

    SetRandom();
    
    int t0=3;
    int l = nBeads/3;

    pimc::levyMove levy(l,0);

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);


    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);
    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal);



    tab.push_back(&close,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&open,0.9,pimc::sector_t::diagonal);

    
    //configurations.setHeadTail(0,M,-1);

    //configurations.setTail(0,-1);
    //configurations.setHead(1,M);
    //configurations.join(0,1);
    //configurations.join(1,0);


    //configurations.fillHeads();


    //Real l2=0;
    //Real l2Var=0;
    //Real l2Error=0;
    int iNewHead=4;

    int nTrials = 10000;
    
    std::ofstream l2Out,x2Out;

    l2Out.open("l2.dat");
    x2Out.open("x2.dat");

    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
        Real l2=0;
        Real l2Var=0;
        Real x2=0;

        accumulate(0,nTrials, [&](
            
            const pimc::configurations_t & confs, const pimc::firstOrderAction & S){

                auto tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;
                //x2+=accumulateX2(0, configurations, 0, M);

            } 
                ,20,pimc::sector_t::offDiagonal
        );


        int nMeasures= nOpen;

        std::cout << nOpen/(nClosed + nOpen) << std::endl;

        l2Out << iBlock << "\t" << l2/nMeasures << std::endl;
        //x2Out << iBlock << "\t" << x2/nMeasures << std::endl;
        
        resetCounters();

    }
    
    l2Out.close();
    x2Out.close();


/* 

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    } */

}


TEST_F(configurationsTest,openChain_harmonic)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpNonInteractingHarmonicAction();
    
    SetRandom();


    SetSeed(time(NULL));

    int t0=10;
    int l = 3;
    int iHead=6;
    int iTail=8;
    

    pimc::levyMove levy(l,0);
    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);


    configurations.setHeadTail(1,M,iTail);
    configurations.setHeadTail(0,iHead,-1);
    configurations.join(1,0);
    configurations.fillHeads();

    pimc::translateMove translate(0.1, 2*M , 0 );


    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.2,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.2,pimc::sector_t::offDiagonal);
    



    int nTrials = 1000;

    std::ofstream x2Out,l2Out;


    l2Out.open("l2.dat");
    x2Out.open("x2.dat");
    
    
    for (int iBlock=0;iBlock < 1000000 ; iBlock++)
    {
            Real l2=0;
            Real x2=0;

            accumulate(
                1000,nTrials, [&](
                
                const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;
                x2+=accumulateX2(0, configurations, 0, M);

                }
                ,100,pimc::sector_t::offDiagonal
                    );
            
            l2Out << iBlock << " " << l2/nOpen << std::endl;
            x2Out << iBlock << " " << x2/nOpen << std::endl;
        
            resetCounters();

    
    }
    
    l2Out.open("l2Out.dat");
    x2Out.open("x2.dat");



/*     l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);

    //EXPECT_NEAR(l2 , 3* (2*M - 1) * timeStep , 2*l2Error);

    std::cout << l2 << " " << l2Error << std::endl;

    for (int d=0;d<getDimensions();d++)
    {
        std::cout << x[d] << " " << x2[d] << std::endl;

    } */

}






TEST_F(configurationsTest,openClosedChain)
{
    Real C=1;
    int nBeads=10;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    SetRandom();


    SetSeed(time(NULL));

    
    int t0=10;
    int l = 3;


    pimc::levyMove levy(l,0);

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);



    tab.push_back(&levy,0.9,pimc::sector_t::diagonal);
    tab.push_back(&levy,0.9,pimc::sector_t::offDiagonal);
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal);
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal);


    tab.push_back(&close,2,pimc::sector_t::offDiagonal);
    tab.push_back(&open,2,pimc::sector_t::diagonal);

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);


    Real l2=0;
    Real l2Var=0;
    Real l2Error=0;

    int nTrials = 1000000;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::diagonal);

    l2/=nClosed;
    l2Var/=nClosed;

    l2Error = std::sqrt((l2Var - l2*l2)/nClosed);

    EXPECT_NEAR(l2 , 3* (M - 1) * timeStep , 3*l2Error);


    resetCounters();


    l2=0;
    l2Var=0;

    accumulate(1000,nTrials, [&](const pimc::configurations_t & confs, const pimc::firstOrderAction & S){Real tmp=accumulateAverageLengthSquare(0,confs) ;l2+=tmp;l2Var+=tmp*tmp;}  ,100,pimc::sector_t::offDiagonal);

    l2/=nOpen;
    l2Var/=nOpen;

    l2Error = std::sqrt((l2Var - l2*l2)/nOpen);


    EXPECT_NEAR(l2 , 3* ( M - l ) * timeStep , 3*l2Error);

    
}




TEST_F(configurationsTest,advanceRecedeGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int l =40;
    int nBeads= 100;

    SetUp(1,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    
    int iHead=60;

    configurations.setHeadTail(0,iHead,-1);

    pimc::advanceHead advanceMove(l,0);
    pimc::recedeHead recedeMove(l,0);

    advanceMove.setFixedLength();
    recedeMove.setFixedLength();

    int nTrails = 100;


    std::array<Real ,getDimensions()> x{DLIST(0,0,0) }; 
    int nMeasurements = 0;
    int i = iHead + l <= M ? iHead + l  : iHead + l - M;
    int iChainNewHead = iHead + l <= M ? 0  : 1;

    const auto & data = configurations.dataTensor();

    SetSeed(56);

    for (int n=0;n<nTrails;n++)
    {
        int currentHead = configurations.getChain(0).head;
        
        if (currentHead == iHead)
        {
            advanceMove.attemptMove(configurations,S,randG);
        }
        else
        {
            for (int d=0;d<getDimensions();d++)
            {
                x[d]+=data(iChainNewHead,d,i);
            }
            nMeasurements+=1;
            recedeMove.attemptMove(configurations,S,randG);
        }

    }

    ASSERT_NEAR( nMeasurements * 1./nTrails, 0.5 , 1e-3 ) ;

    for (int d=0;d<getDimensions();d++)
    {
        x[d]/=nMeasurements;
        ASSERT_NEAR ( data(0,d,iHead)  , x[d] , 1e-2 );

    }



}

TEST_F(configurationsTest,swapGrandCanonical_distributionReconstructedChain)
{
    Real C=1;
    int l =40;
    int nBeads= 100;

    SetUp(2,nBeads,1);
    SetGrandCanonicalEnsamble(0);
    SetUpFreeParticleAction();
    

    int iHead=30;

    configurations.setHeadTail(0,iHead,-1);

    pimc::swapMove swap(l,2,0);

    swap.setFixedLength();

    SetSeed(58);

    int nTrails = 1000000;
    int nBurns = 1000;
    int nMeasurements=0;

    const auto & data = configurations.dataTensor();

    pimc::configurations_t configurationsInitial(configurations);
    const auto & dataInitial = configurations.dataTensor();

    

    std::array<Real ,getDimensions()> x{DLIST(0,0,0) }; 
    int i = iHead + 10;

    auto meanExpected=meanBeadFixedLengths(0 ,  1 , iHead, iHead + l, i , timeStep,configurations);


    for (int n=0;n<nBurns;n++)
    {
        bool accept=swap.attemptMove(configurations,S,randG);
    }
    for (int n=0;n<nTrails;n++)
    {
        bool accept=swap.attemptMove(configurations,S,randG);

        assert( configurations.getChain(0).hasHead() or configurations.getChain(1).hasHead( ));


        if ( not configurations.getChain(0).hasHead()  )
        {
            
            for (int d=0;d<getDimensions();d++)
            {
                x[d]+=data(0,d,i);
            }
            nMeasurements+=1;
        }

    }

    for (int d=0;d<getDimensions();d++)
            {
                x[d]/=nMeasurements;
            }
    
    

    for (int d=0;d<getDimensions();d++)
    {

       

        int iChainHead = configurations.getGroups()[0].heads[0];
        int iChainCont = iChainHead == 0 ? 1 : 0;

            for (int t=0;t<=iHead;t++)
            {
                assert(data(0,d,t) == dataInitial(0,d,t)    );
                assert(data(1,d,t) == dataInitial(1,d,t)    );
            }

            for (int t=iHead+l;t<=M;t++)
            {
                assert(data(iChainCont,d,t) == dataInitial(1,d,t)    );
            }
         ASSERT_NEAR(x[d],meanExpected[d],1e-2);
         
           }

}


TEST(moves,openCloseGrandCanonical)
{   
    int N=1;
    int M=100;
    Real Beta = 1;

    int seed=812; // 356
    int buffer=2;

    int nChains=N + buffer;

    Real timeStep=Beta/M;
    pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});


    configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations.setChemicalPotential({0});

    pimc::geometryPBC_PIMC geo(300,300,300);

    
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    std::shared_ptr<pimc::action> sV= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);


    randomGenerator_t randG(seed);

   

    pimc::firstOrderAction S(sT,  sV);

    auto & data=configurations.dataTensor();

    for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=uniformDistribution(randG);
            }
    
    Real C = 1;
    int l=30;
    int t0=6;
    int t1=t0 + l;
    int i = t0 + 20;



    pimc::tableMoves tab;

    pimc::openMove open(C, 0, l );
    pimc::closeMove close(C, 0, l );

    open.setStartingBead(t0);
    close.setStartingBead(t0);

    open.setLengthCut(l);
    close.setLengthCut(l);

    tab.push_back(& open,1,pimc::sector_t::diagonal);
    tab.push_back(& close, 1, pimc::sector_t::offDiagonal);


    configurations.fillHeads();

    Real nTrials = 100000;
    Real nBurns = 100;

    Real nOpen=0;

    std::array<Real,getDimensions() > x{DLIST(0,0,0)};
    std::array<Real ,getDimensions() > x2{DLIST(0,0,0)};


    for(int k=0;k<nBurns;k++)
    {
        bool accept= tab.attemptMove(configurations,S,randG);    
    }
    
    for(int k=0;k<nTrials;k++)
    {

        if ( configurations.isOpen() )
        {
            nOpen+=1;
        }
        else
        {
                for(int d=0;d<getDimensions();d++)
                {
                    x[d]+=data( 0 ,d, i) ;
                    x2[d]+=std::pow(data(0,d,i),2);
                }
        }

        bool accept= tab.attemptMove(configurations,S,randG);

    }

    Real nClosed = nTrials - nOpen;

    Real openRatio = nOpen/nTrials;
    
    std::array<Real,getDimensions()> meanExpected;
    std::array<Real,getDimensions()> varianceExpected;

    Real D = 0.5;
    Real mass = 1;

    std::array<Real,getDimensions()> difference;

    for (int d=0;d<getDimensions();d++)
    {
        difference[d]= data(0,d,t0) - data(0,d,t1);

        varianceExpected[d]=1./(1./(i - t0) + 1./(t1 - i) )* 2 * D * timeStep / mass;

        meanExpected[d]=(data(0,d,t0)/ (i-t0) + 
        data(0,d,t1)/( t1 - i) )/(1./(i - t0) + 1./(t1 - i) );

        x[d]/=nClosed;
        x2[d]/=nClosed;


    }


    Real expectedOpenRatio= 1/( 1 + exp(pimc::freeParticleLogProbability(difference,S.getTimeStep()*l,mass))/C );

    std::cout << "Open ratio error : " << std::abs(openRatio - expectedOpenRatio)/expectedOpenRatio <<  std::endl;

    ASSERT_NEAR((openRatio - expectedOpenRatio)/expectedOpenRatio , 0 , 1e-2);


    std::cout << "Mean error: " << std::abs((x[0]  -  meanExpected[0])/meanExpected[0])<< std::endl;
    ASSERT_NEAR( (x[0]  -  meanExpected[0])/meanExpected[0] , 0 , 1e-2);

    std::cout << "var error : " << std::abs(x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] << std::endl;

    ASSERT_NEAR(( x2[0] - x[0]*x[0] - varianceExpected[0] )/varianceExpected[0] , 0 , 1e-2);

}


TEST(run,free_harmonic_oscillator_grandCanonical)
{   
    int N=1;
    int M=20;
    Real Beta = 1;


    std::cout << "dimensions" << getDimensions() << std::endl;
    int seed= time(NULL);
    //int seed= 42;

    int buffer=1000;


    int nChains=N + buffer;

    std::srand((unsigned int) seed);
    randomGenerator_t randG(seed);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);


    pimc::geometryPBC_PIMC geo(3000,3000,3000);


    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    Real mu = getDimensions()/2. * 0.0;

    //Real mu =  0 ;


    configurations.setEnsamble(pimc::ensamble_t::grandCanonical);
    configurations.setChemicalPotential( {mu} );

    ASSERT_DOUBLE_EQ(configurations.getChemicalPotential(0) , mu );

    auto & data=configurations.dataTensor();

    for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=uniformDistribution(randG);
            }


    //std::cout << configurations.dataTensor().sum();
    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();

    pimc::levyReconstructor reconstructor(M);

    

    Real delta=0.1;

    pimc::translateMove translMove(delta,(M+2)*(nChains),0);

    Real C = 1e-4;
    int l =16;
    int lShort=10;

    int t0=8;

    pimc::levyMove freeMoves(l,0);

    pimc::openMove openMove(C,0,lShort);
    pimc::closeMove closeMove(C,0,lShort);


    //openMove.setStartingBead(t0);
    //closeMove.setStartingBead(t0);

    //openMove.setLengthCut(4);
    //closeMove.setLengthCut(4);

    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);


    pimc::advanceHead advanceHeadMove(lShort,0);
    pimc::recedeHead recedeHeadMove(lShort,0);

    pimc::advanceTail advanceTailMove(l,0);
    pimc::recedeTail recedeTailMove(l,0);

    pimc::createWorm addMove(C,0,lShort,1);
    pimc::deleteWorm removeMove(C,0,lShort,1);


    pimc::swapMove swapMove(lShort,N,0);


    pimc::tableMoves table;

    table.push_back(& freeMoves,0.7,pimc::sector_t::diagonal,"levy");
    table.push_back(& openMove,0.1,pimc::sector_t::diagonal,"open");
    table.push_back(& addMove,0.1,pimc::sector_t::diagonal,"addMove");
    table.push_back(& translMove,0.1,pimc::sector_t::diagonal,"translate");
    
    
    
    table.push_back(& freeMoves,0.2,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& translMove,0.1,pimc::sector_t::offDiagonal,"translate");
    table.push_back(& closeMove,0.1,pimc::sector_t::offDiagonal,"close");
    table.push_back(& removeMove,0.1,pimc::sector_t::offDiagonal,"removeMove");
    table.push_back(& moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    table.push_back(& advanceHeadMove,0.1,pimc::sector_t::offDiagonal,"advanceHead");
    table.push_back(& recedeHeadMove,0.1,pimc::sector_t::offDiagonal,"recedeHead");
    //table.push_back(& advanceTailMove,0.1,pimc::sector_t::offDiagonal,"advanceTail");
    //table.push_back(& recedeTailMove,0.1,pimc::sector_t::offDiagonal,"recedeTail");
    table.push_back(& swapMove,0.1,pimc::sector_t::offDiagonal,"swap");



    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );


    Real R0=0.1;
    Real V0=1;



  
      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-(r*r));} ,
         [=](Real r) {return -2*r*V0*exp(-r*r)  ;}
          );




    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,nChains,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};
    
    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    std::shared_ptr<pimc::action> sEmpty= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);

    pimc::firstOrderAction S(sT,  sV);

    
    int nTimes = 100000;
    int success = 0;
    int subSteps=10000;
    int correlationSteps=10;
   
    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(nChains, M);

    pimc::particleNumberEstimator particleEstimator(0);

    Real e=0;
    Real e2=0;
    Real sumN=0;


    std::ofstream f;
    std::ofstream fV;
    std::ofstream fN;
    

    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    configurations.fillHeads();

    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat");
    fV.open("energyVirial.dat");
    fN.open("N.dat");

    int iSuspicious=0;
    int iOpen=0;
    int nBurns=0 * subSteps * correlationSteps;

    for (int i=0;i< nBurns ; i++)
    {
        bool accepted=table.attemptMove(configurations, S, randG);

    }
            
    


    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0,nStep=0;
        int nMeasurements=0;

        for (int k=0;k< subSteps;k++)
        {
            
            for (int j=0;j<correlationSteps;j++)
            {
                bool accepted=table.attemptMove(configurations, S, randG);

                if (accepted)
                {success+=1;}
            }
            
            if (!configurations.isOpen() )
            {
                Real tmp=energyEstimator(configurations,S);
                Real tmp1=viriralEnergy(configurations,S);

                Real tmp2=particleEstimator(configurations,S);
            
                nMeasurements++;

                /* if (tmp < -6)
                {
                    std::string name="configurations/sampleSusp"+std::to_string(iSuspicious);
                    configurations.save(name,"pdb"); 
                    iSuspicious++;             
                } */
                eStep+=tmp;
                eVirialStep+=tmp1;
                nStep+=tmp2;
            }
            else
            {
                //configurations.save("configurations/sampleOpen"+std::to_string(iOpen),"pdb");
                //iOpen++;
            }
            
        }

        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        fN << i + 1 << "\t" << nStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;
        std::cout << "N: " << configurations.nParticles() << std::endl;

        table >> std::cout;

        //configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
    }

    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}





TEST(run,free)
{   

    Real density=0.15884256651199277;
    int N=10;

    int M=10;
    Real Beta = 1;
    Real lBox = std::pow(N/density,1./3) ;
    

    pimc::geometryPBC_PIMC geo(lBox,lBox,lBox);


    Real timeStep = Beta/M;

    pimc::particleGroup groupA{ 0 , N-1, N - 1 , 1.0};
    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor(); 
    data.setRandom();

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<M;j++)
        for(int d=0;d<getDimensions();d++)
        {
            data(i,d,j)=(data(i,d,j)-0.5 )*lBox;
        }
    }
    


    //configurations.join(0,1);
    //configurations.join(1,0);
    
    configurations.fillHeads();

    //configurations.dataTensor()(0,0,M)=configurations.dataTensor()(0,0,0) + 5*lBox;
    


    pimc::levyReconstructor reconstructor(M);

    pimc::levyMove freeMoves(M/3,0);


    Real delta=0.1;

    pimc::translateMove translMove(delta,(M+1)*N,0);

    Real C = 1e-1;
    int l = 4;


    pimc::openMove openMove(C,0,l);
    pimc::closeMove closeMove(C,0,l);

    pimc::moveHead moveHeadMove(l,0);
    pimc::moveTail moveTailMove(l,0);

    pimc::swapMove swapMove(4,N,0);

    pimc::tableMoves table;

    table.push_back(& freeMoves,0.8,pimc::sector_t::offDiagonal,"levy");
    table.push_back(& freeMoves,0.8,pimc::sector_t::diagonal,"levy");

    table.push_back(& translMove,0.2,pimc::sector_t::diagonal,"translate");
    table.push_back(& translMove,0.2,pimc::sector_t::offDiagonal,"translate");

    table.push_back(& openMove,0.1,pimc::sector_t::diagonal,"open");
    table.push_back(& closeMove,0.1,pimc::sector_t::offDiagonal,"close");

    table.push_back(& moveHeadMove,0.4,pimc::sector_t::offDiagonal,"moveHead");
    table.push_back(& moveTailMove,0.4,pimc::sector_t::offDiagonal,"moveTail");

    //table.push_back(& swapMove,1.9,pimc::sector_t::offDiagonal,"swap");


    randomGenerator_t randG(368);

     std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


   
     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0 ;} ,
         [](Real r) {return 0  ;}
         );

    Real alpha=1/std::pow (0.1*lBox,2);

    Real V0=100;


      auto V2 = pimc::makeIsotropicPotentialFunctor(
         [=](Real r) {return V0*exp(-alpha*(r*r));} ,
         [=](Real r) {return -2*r*V0*alpha*exp(-alpha*r*r)  ;}
          );


    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    std::shared_ptr<pimc::action>  sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V2)>  >(timeStep,N,M,V2 ,geo,0,0);    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,sV2B};
    
    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);


    /*
    std::shared_ptr<pimc::action> sV=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo);
    */

    pimc::firstOrderAction S(sT,  sV);
    int nTimes = 1000;
    int success = 0;
    int subSteps=1000;
    int correlationSteps=100;

    pimc::thermodynamicEnergyEstimator energyEstimator;

    pimc::virialEnergyEstimator viriralEnergy(N, M);


    Real e=0;
    Real e2=0;

    std::ofstream f;
    std::ofstream fV;


    if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }


    configurations.save("configurations/sample"+std::to_string(0),"pdb");

    f.open("energy.dat");
    fV.open("energyVirial.dat");

    for (int i=0;i< nTimes ; i++)
    {
        Real eStep=0,eVirialStep=0;
        int nMeasurements=0;

        for (int k=0;k< subSteps;k++)
        {
            
            for (int j=0;j<correlationSteps;j++)
            {
                bool accepted=table.attemptMove(configurations, S, randG);

                if (accepted)
                {success+=1;}
            }
            
            if (!configurations.isOpen() )
            {
                Real tmp=energyEstimator(configurations,S);
                Real tmp1=viriralEnergy(configurations,S);
            
                nMeasurements++;
            
                eStep+=tmp;
                eVirialStep+=tmp1;
            }
            else
            {
               
            }
            
        }


        f << i + 1 << "\t" << eStep/nMeasurements << std::endl ;
        fV << i + 1 << "\t" << eVirialStep/nMeasurements << std::endl ;

        //std::cout << e << std::endl;
        std::cout << "Energy: " << eStep/nMeasurements << std::endl;
        std::cout << "Acceptance ratio: " << success*1./((i+1)*subSteps*correlationSteps) << std::endl;

        table >> std::cout;

        configurations.save("configurations/sample"+std::to_string(i+1),"pdb");
    }
    
    f.close();
    ASSERT_TRUE( (success*1./nTimes )> 0);
    std::cout << "END." << std::endl;
}





















TEST(configurations, forces)
{
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N + 1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();
    const auto & tags = configurations.getTags();

    data.setRandom();

    configurations.fillHeads();

    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );

    pimc::totalForce totalForce(M,N);


    auto oneBodyForce = 
    std::make_shared<pimc::oneBodyIsotropicForce<decltype(V) > >( V , 0) ;


    totalForce.add(oneBodyForce);

    totalForce.evaluate(configurations,{0,M-1},{0,N-1});


    const auto & forceData = totalForce.getForceData();


    for (int t=0;t<M;t++)
    {
        for (int i=0;i<N;i++)
        {
            Real r2=0;
            for(int d=0;d<getDimensions();d++)
            {
                r2+=data(i,d,t)*data(i,d,t);
            }    
            Real r = std::sqrt(r2);

            for(int d=0;d<getDimensions();d++)
            {
            
                ASSERT_NEAR( forceData(i,d,t)   , V.radialDerivative(r)*data(i,d,t)/r  , 1e-6); 
            }
        }
    }



}

#include "../pimc/potentials.h"


TEST(configurations, potentials_canonical)
{
    const int N = 10;
    const int M = 50;

    pimc::particleGroup groupA{ 0 , N-1, N + 1 , 1.0};

    pimc::pimcConfigurations configurations(M , getDimensions() , {groupA});

    auto & data = configurations.dataTensor();
    const auto & tags = configurations.getTags();

    data.setRandom();

    configurations.fillHeads();

    auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );

    pimc::potentials totalPotential;

    auto oneBodyPotential = 
    std::make_shared<pimc::oneBodyIsotropicPotentialSecondOrder<decltype(V) > >( V , 0) ;

    totalPotential.add(oneBodyPotential);

    Real sum=totalPotential.evaluate(configurations,{0,M-1},{0,N-1});

    Real sumCheck=0;

    for (int t=0;t<M;t++)
    {
        Real prefactor = ( ( t== 0) or (t==M) ) ? 0.5 : 1 ;
        for (int i=0; i<N ;i++)
        {
            Real r2=0;
            for( int d=0 ; d<getDimensions(); d++)
            {
                r2+=data(i,d,t)*data(i,d,t);
            }    
            Real r = std::sqrt(r2);
            sumCheck+=prefactor*V(r);
        }
    }


    ASSERT_NEAR(  sumCheck , sum , 1e-6 );


    // chin action potential

    std::array<Real,2> coeff = {1,2};


    auto oneBodyPotentialChin = 
    std::make_shared<pimc::oneBodyIsotropicPotentialChin<decltype(V) > >( V , 0, coeff) ;




}



