#include "integrationTests.h"

#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"
#include "testConfigurations.h"
#include "../pimc/pimcObservables.h"
#include <filesystem>

namespace fs = std::filesystem;


TEST_F( configurationsTest, semiOpenMove)
{

    int nBeads=10;
    int N=10;
    Real beta=0.1* nBeads;
    int nTrials=1000;

    std::array<Real,3> lBox{ TRUNCATE_D(1,1,1)};
    SetUp(N,nBeads,beta , lBox );
    SetUpFreeParticleAction();

    SetGrandCanonicalEnsamble( 0 );

    SetSeed( 145) ;
    SetSeed(456);

    SetRandom(lBox);

    configurations.fillHeads();
    
    int t0=3;
    int l = 6;
    int iChain=6;

    int nBurnsIn=10000;



    auto j = R"(
        {
            "C" : 1e-2,
            "setA": 0,
            "setB" : 1,
            "reconstructionMaxLength": 6
        }
            )"_json;
    
    pimc::semiOpenMove openMove(j);

    pimc::levyMove levyMove(l,0);
    
    std::cout << "Equilibrating.." << std::endl;
    for(int n=0;n<nBurnsIn;n++)
    {
        levyMove.attemptMove(configurations,S,randG);
    }

    openMove.setStartingBead(t0);
    openMove.setStartingChain(iChain);

    const auto & data = configurations.dataTensor();

    
    std::cout << "Test openings.." << std::endl;

    Real nAccepted=0;
    for(int t=0;t<nTrials;t++)
    {
        pimc::configurations_t configurations2(configurations);
        const auto & data2 = configurations.dataTensor();
        
        bool accept=openMove.attemptMove( configurations2,S,randG);
        if (accept) nAccepted++;

        if (accept)
        {
            ASSERT_TRUE( configurations2.isOpen(0) );
            TESTCHAIN(configurations2,iChain,M,t0-1);

            int iChainNext=configurations2.getChain(iChain).next;
            ASSERT_NE(iChainNext,-1);
            TESTCHAIN(configurations2,iChainNext,t0,-1);
        }
        else
        {
            ASSERT_FALSE( configurations2.isOpen(0) );
            TESTCHAIN(configurations2,iChain,M,-1);

            for(int t=0;t<M;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    ASSERT_NEAR( data(iChain,d,t),data2(iChain,d,t) ,1e-9);
                }

            }


        }
    }



    std::cout << "Acceptance ratio: " <<  nAccepted * 1./nTrials << std::endl;

    
    ASSERT_GT(nAccepted,0);
}

TEST_F( configurationsTest, semiOpenCloseMove )
{

    auto j= R"(
        {
        "C": 1e-1,
        "l" : 4,
        "t0" : 8,
        "nBlocks" : 10000,
        "stepsPerBlock" : 100000
        }
    )"_json;


    if (arguments::testInputFile != "" )
    {
        std::ifstream ifs(arguments::testInputFile);
        j = pimc::json_t::parse(ifs);
    }



    Real C=j["C"].get<Real>();
    int l=j["l"].get<int>();
    int t0=j["t0"].get<int>();
     int nBlocks=j["nBlocks"].get<int>();
    int nStepsPerBlock=j["stepsPerBlock"].get<int>();

    auto jOpenClose = R"(
        {
            "C" : 1e-1,
            "setA": 0,
            "setB" : 1,
            "reconstructionMaxLength": 4
        }
            )"_json;
    
    
    jOpenClose["C"]=C;
    jOpenClose["reconstructionMaxLength"]=l;

    std::cout << arguments::testInputFile << std::endl;

    int nBeads=10;
    int N=1;
    Real beta=0.1* nBeads;
   

    std::array<Real,3> lBox{ TRUNCATE_D(1,1,1)};
    SetUp(N,nBeads,beta , {30000,300000,300000} );
    SetUpNonInteractingHarmonicAction();
    SetGrandCanonicalEnsamble( 0 );

   // SetSeed(456);

    SetRandom(lBox);

    configurations.fillHeads();
    int lLong=int(M*0.6);
    int lShort=int(M*0.3);

    int iChain=0;    


    int nBurnsIn=10000;
    
    pimc::semiOpenMove openMove(jOpenClose);
    pimc::semiCloseMove closeMove(jOpenClose);
    pimc::levyMove levyMove(lLong,0);

    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveHead moveTailMove(lShort,0);
    
    Real delta=0.05;
    pimc::translateMove translMove(delta,(M+1)*(N+1),0);



    tab.push_back(&openMove,0.2,pimc::sector_t::diagonal,"open");
    tab.push_back(&levyMove,0.7,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translMove,0.1,pimc::sector_t::diagonal,"translate");

    tab.push_back(&closeMove,0.2,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&levyMove,0.6,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&moveHeadMove,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMove,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translMove,0.1,pimc::sector_t::offDiagonal,"translate");

    //configurations.setHeadTail(0,  M , t0-1);
    //configurations.setHeadTail(1,  t0 , -1);
    //configurations.join(0,1);

    configurations.fillHeads();


   auto lengthEstimator=std::make_shared< pimc::lengthEstimator >( M  );

   lengthEstimator->setStartingBead(0,t0);


    auto lengthObClosed = std::make_shared<pimc::scalarObservable>(lengthEstimator,std::string("lengthClosed"),false);

    auto lengthObOpen = std::make_shared<pimc::scalarObservable>(lengthEstimator,std::string("lengthOpen"),false);


    std::cout << "Equilibrating.." << std::endl;


    for(int n=0;n<nBurnsIn;n++)
    {
        levyMove.attemptMove(configurations,S,randG);
    }

    openMove.setStartingBead(t0);
    openMove.setStartingChain(iChain);
    openMove.setLength(l);

    closeMove.setStartingBead(t0);
    closeMove.setStartingChain(iChain);
    closeMove.setLength(l);

    const auto & data = configurations.dataTensor();

    std::array<Real,getDimensions()> length{0,0,0};

    Real nOpen=0;
    Real nClosed=0;
    Real ratio=0;
    int n=0;
    
    std::ofstream openRatioOut("openRatio.dat");

    std::cout << "Start sampling..." << std::endl;
    for(int t=0;t<nBlocks;t++)
    {
       while ( (nOpen < nStepsPerBlock) and (nClosed < nStepsPerBlock) and ( n < nStepsPerBlock ) )
        {
            n+=1;
            bool accept=tab.attemptMove(configurations,S,randG);
            if (configurations.isOpen() )
            {
                nOpen++;
                ratio+=1;

                lengthObOpen->accumulate(configurations,S);

            }
            else
            {
                lengthObClosed->accumulate(configurations,S);
                nClosed++;


            }
        }
    


    std::cout << "Open sector ratio: "<< nOpen*1/(nOpen + nClosed) << std::endl;
     tab >> std::cout ;

    if (nOpen == nStepsPerBlock)
        {
            lengthObOpen->out(t);
            lengthObOpen->clear();
            nOpen=0;
        }
    
    if ( nClosed== nStepsPerBlock)
        {
            lengthObClosed->out(t);
            lengthObClosed->clear();
            nClosed=0;

        }

    if ( n == nStepsPerBlock)
        {
            ratio/=1.*n;
            openRatioOut << t << " " << ratio << std::endl;
            ratio=0;
            n=0;
        }

    }

}


TEST_F( configurationsTest, advanceRecedeHeadTail_run )
{
    auto j= R"(
        {
        "C": 1e-1,
        "l" : 2,
        "t0" : 6,
        "nBlocks" : 10000,
        "stepsPerBlock" : 100000
        }
    )"_json;


    if (arguments::testInputFile != "" )
    {
        std::ifstream ifs(arguments::testInputFile);
        j = pimc::json_t::parse(ifs);
    }

    Real C=j["C"].get<Real>();
    int l=j["l"].get<int>();
    int t0=j["t0"].get<int>();
     int nBlocks=j["nBlocks"].get<int>();
    int nStepsPerBlock=j["stepsPerBlock"].get<int>();


    auto jAdvanceHeadTail = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4
        }
            )"_json;
    
    jAdvanceHeadTail["reconstructionMaxLength"]=l;

    std::cout << arguments::testInputFile << std::endl;


    int nBeads=10;
    int NA=2;
    int NB=2;

    Real beta=0.1* nBeads;
    std::array<Real,3> lBox{ TRUNCATE_D(1,1,1)};
    SetUp( {NA,NB}, nBeads,beta , {30000,300000,300000} );
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetGrandCanonicalEnsamble( {0,0} );

    const auto & groupA=configurations.getGroups()[0];
    const auto & groupB=configurations.getGroups()[1];


   // SetSeed(456);

    SetRandom(lBox);

    configurations.fillHeads();
    
    int lLong=int(M*0.6);
    int lShort=int(M*0.3);

    int iChain=0;    


    int nBurnsIn=10000;
    
    pimc::advanceHeadTail advanceMove(jAdvanceHeadTail);
    pimc::recedeHeadTail recedeMove(jAdvanceHeadTail);

    advanceMove.setFixedLength();
    recedeMove.setFixedLength();


    pimc::levyMove levyMoveA(lLong,0);
    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveHead moveTailMoveA(lShort,0);

    pimc::levyMove levyMoveB(lLong,1);
    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveHead moveTailMoveB(lShort,1);

    Real delta=0.5;
    pimc::translateMove translMoveA(delta,(M+1)*configurations.nChains()*2,0);
    pimc::translateMove translMoveB(delta,(M+1)*configurations.nChains()*2,1);


    tab.push_back(&levyMoveA,0.8,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::offDiagonal,"translateA");

    tab.push_back(&levyMoveB,0.8,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.05,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::offDiagonal,"translateB");


    tab.push_back(&levyMoveA,0.9,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::diagonal,"translateA");

    tab.push_back(&levyMoveB,0.9,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::diagonal,"translateB");

    //configurations.setHeadTail(0 ,  M-1 , t0-1);

    int deltaL=0;
    configurations.setHeadTail(0 + groupA.iStart,  M , t0 + deltaL -1);
    configurations.setHeadTail(1 + groupA.iStart,  t0  , -1);
    configurations.join(0 + groupA.iStart,1 + groupA.iStart);

    configurations.setHeadTail(0 + groupB.iStart,  M , t0 -1);

    configurations.setHeadTail(1 + groupB.iStart,  t0 + deltaL , -1);
    configurations.join(0 + groupB.iStart,1 + groupB.iStart);

    //configurations.setHeadTail(1,  t0 , -1);
    //configurations.join(0,1); */

    configurations.fillHeads();


   auto lengthEstimatorAShort=std::make_shared< pimc::lengthEstimator >( );
   auto lengthEstimatorALong=std::make_shared< pimc::lengthEstimator >( );
   auto lengthEstimatorBShort=std::make_shared< pimc::lengthEstimator >();
   auto lengthEstimatorBLong=std::make_shared< pimc::lengthEstimator >( );

    

    lengthEstimatorAShort->setStartFromHeadOrTail(0);
    lengthEstimatorBShort->setStartFromHeadOrTail(1);
    lengthEstimatorALong->setStartFromHeadOrTail(0);
    lengthEstimatorBLong->setStartFromHeadOrTail(1);


    auto lengthAObShort = std::make_shared<pimc::scalarObservable>(lengthEstimatorAShort,std::string("lengthAShort"),false);
    auto lengthAObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorALong,std::string("lengthALong"),false);

    auto lengthBObShort = std::make_shared<pimc::scalarObservable>(lengthEstimatorBShort,std::string("lengthBShort"),false);
    auto lengthBObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorBLong,std::string("lengthBLong"),false);

    
    std::cout << "Equilibrating..." << std::endl;

    for(int n=0;n<nBurnsIn;n++)
    {
        tab.attemptMove(configurations,S,randG);
    }

    const auto & data = configurations.dataTensor();

    Real nShort=0;
    Real nLong=0;
    Real ratio=0;
    Real n=0;

    std::ofstream openRatioOut("openRatio.dat");

    std::cout << "Start sampling..." << std::endl;
    for(int t=0;t<nBlocks;t++)
    {
       while ( (nShort < nStepsPerBlock) and (nLong < nStepsPerBlock) and ( n < nStepsPerBlock ) )
        {
            n+=1;
            bool accept=tab.attemptMove(configurations,S,randG);

            if ( configurations.getChain(groupA.iStart + 1).head == t0 )
            {
                advanceMove.attemptMove(configurations,S,randG);
            }
            else 
            {
                recedeMove.attemptMove(configurations,S,randG);
            }

            if (  configurations.getChain(groupA.iStart + 1).head == t0 )
            {
                nShort++;
                ratio+=1;
                lengthAObShort->accumulate(configurations,S);
                lengthBObShort->accumulate(configurations,S);
            }
            else
            {
                lengthAObLong->accumulate(configurations,S);
                lengthBObLong->accumulate(configurations,S);
                nLong++;
            }
        }

    std::cout << "Short length fraction: "<< nShort*1./n << std::endl;

    tab >> std::cout ;

    if (nShort == nStepsPerBlock)
        {
            lengthBObShort->out(t);
            lengthBObShort->clear();

            lengthAObShort->out(t);
            lengthAObShort->clear();

            nShort=0;
        }
    
    if ( nLong== nStepsPerBlock)
        {
            lengthBObLong->out(t);
            lengthBObLong->clear();

            lengthAObLong->out(t);
            lengthAObLong->clear();

            nLong=0;
        }

    if ( n == nStepsPerBlock)
        {
            ratio/=1.*n;
            openRatioOut << t << " " << ratio << std::endl;
            ratio=0;
            n=0;
        }


    }

}


TEST_F( configurationsTest, semiFullOpenClose_run )
{

    auto j= R"(
        {
        "CA": 1e-1,
        "l" : 4,
        "t0" : 6,
        "nBlocks" : 10000,
        "stepsPerBlock" : 100000
        }
    )"_json;


    if (arguments::testInputFile != "" )
    {
        std::ifstream ifs(arguments::testInputFile);
        j = pimc::json_t::parse(ifs);
    }

    Real CA=j["CA"].get<Real>();
    int l=j["l"].get<int>();
    int t0=j["t0"].get<int>();
     int nBlocks=j["nBlocks"].get<int>();
    int nStepsPerBlock=j["stepsPerBlock"].get<int>();

    auto jFullSemiCanonicalOpenSettings = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4,
            "C" : 1
        }
            )"_json;
    
    jFullSemiCanonicalOpenSettings["reconstructionMaxLength"]=l;
    jFullSemiCanonicalOpenSettings["C"]=CA;


    int nBeads=10;
    int NA=1;
    int NB=2;


    Real beta=0.1* nBeads;
    std::array<Real,3> lBox{ TRUNCATE_D(1,1,1)};
    SetUp( {NA,NB}, nBeads,beta , {30000,300000,300000} );
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetGrandCanonicalEnsamble( {0,0});


    const auto & groupA=configurations.getGroups()[0];
    const auto & groupB=configurations.getGroups()[1];


   // SetSeed(456);

    SetRandom(lBox);

    configurations.fillHeads();
    
    int lLong=int(M*0.6);
    int lShort=int(M*0.3);

    int iChain=0;    


    int nBurnsIn=10000;
    
    pimc::fullSemiCanonicalOpenMove openMove(jFullSemiCanonicalOpenSettings);
    pimc::fullSemiCanonicalCloseMove closeMove(jFullSemiCanonicalOpenSettings);

    openMove.setLength(l);
    closeMove.setLength(l);




    pimc::levyMove levyMoveA(lLong,0);
    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveHead moveTailMoveA(lShort,0);

    pimc::levyMove levyMoveB(lLong,1);
    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveHead moveTailMoveB(lShort,1);

    Real delta=0.5;
    pimc::translateMove translMoveA(delta,(M+1)*configurations.nChains()*2,0);
    pimc::translateMove translMoveB(delta,(M+1)*configurations.nChains()*2,1);    

    tab.push_back(&levyMoveA,0.8,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::offDiagonal,"translateA");

    tab.push_back(&levyMoveB,0.8,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.05,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::offDiagonal,"translateB");


    tab.push_back(&levyMoveA,0.9,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::diagonal,"translateA");

    tab.push_back(&levyMoveB,0.9,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::diagonal,"translateB");

    //configurations.setHeadTail(0 ,  M-1 , t0-1);

    configurations.setHeadTail(0 + groupB.iStart,  M , t0-1);
    configurations.setHeadTail(1 + groupB.iStart,  t0  , -1);
    configurations.join(0 + groupB.iStart,1 + groupB.iStart);


    /* configurations.setHeadTail(0 + groupA.iStart,  t0 , t0 + l - M -1);

    configurations.setHeadTail(0 + groupB.iStart,  M , t0 -1);
    configurations.setHeadTail(2 + groupB.iStart,  t0 + l - M , - 1);
    
    configurations.join(groupB.iStart + 0 ,groupB.iStart + 1);
    configurations.join(groupB.iStart + 1 ,groupB.iStart + 2); */

    /* configurations.setHeadTail(0 + groupB.iStart,  M , t0-1);
    configurations.setHeadTail(1 + groupB.iStart,  t0 + l  , -1);
    configurations.join(0 + groupB.iStart,1 + groupB.iStart);

    configurations.setHeadTail(0 + groupA.iStart,  M , t0 + l-1);
    configurations.setHeadTail(1 + groupA.iStart,  t0  , -1);
    configurations.join(0 + groupA.iStart,1 + groupA.iStart); */


   configurations.fillHeads();

   
   auto lengthEstimatorBShort=std::make_shared< pimc::lengthEstimator >( M );
   auto lengthEstimatorBLong=std::make_shared< pimc::lengthEstimator >( M + l );
   auto lengthEstimatorAShort=std::make_shared< pimc::lengthEstimator >( M );
   auto lengthEstimatorALong=std::make_shared< pimc::lengthEstimator >( M - l );

    lengthEstimatorBShort->setStartFromHeadOrTail(1);
    lengthEstimatorAShort->setStartingBead(groupA.iStart + 0,t0);
    lengthEstimatorBLong->setStartFromHeadOrTail(1);
    lengthEstimatorALong->setStartFromHeadOrTail(0);
   


    auto lengthAObShort = std::make_shared<pimc::scalarObservable>(lengthEstimatorAShort,std::string("lengthAShort"),false);
    auto lengthAObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorALong,std::string("lengthALong"),false);

    auto lengthBObShort = std::make_shared<pimc::scalarObservable>(lengthEstimatorBShort,std::string("lengthBShort"),false);
    auto lengthBObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorBLong,std::string("lengthBLong"),false);

    std::cout << "Equilibrating..." << std::endl;

    for(int n=0;n<nBurnsIn;n++)
    {
        tab.attemptMove(configurations,S,randG);
    }

    const auto & data = configurations.dataTensor();

    Real nShort=0;
    Real nLong=0;
    Real ratio=0;
    Real n=0;

    std::ofstream openRatioOut("openRatio.dat");

    std::cout << "Start sampling..." << std::endl;
    for(int t=0;t<nBlocks;t++)
    {
       while ( (nShort < nStepsPerBlock) and (nLong < nStepsPerBlock) and ( n < nStepsPerBlock ) )
        {
            n+=1;
            bool accept=tab.attemptMove(configurations,S,randG);

            if ( not configurations.isOpen(0) )
            {
                openMove.attemptMove(configurations,S,randG);
            }
            else 
            {
                closeMove.attemptMove(configurations,S,randG);
            }


            if ( not configurations.isOpen(0) )
            {
                nShort++;
                ratio+=1;
                lengthAObShort->accumulate(configurations,S);
                lengthBObShort->accumulate(configurations,S);
            }
            else
            {
                lengthAObLong->accumulate(configurations,S);
                lengthBObLong->accumulate(configurations,S);
                nLong++;
            }
        }

    std::cout << "Short length fraction: "<< nShort*1./n << std::endl;

    tab >> std::cout ;

    if (nShort == nStepsPerBlock)
        {
            lengthBObShort->out(t);
            lengthBObShort->clear();

            lengthAObShort->out(t);
            lengthAObShort->clear();

            nShort=0;
        }
    
    if ( nLong== nStepsPerBlock)
        {
            lengthBObLong->out(t);
            lengthBObLong->clear();

            lengthAObLong->out(t);
            lengthAObLong->clear();

            nLong=0;
        }

    if ( n == nStepsPerBlock)
        {
            ratio/=1.*n;
            openRatioOut << t << " " << ratio << std::endl;
            ratio=0;
            n=0;
        }


    }

}

TEST_F( configurationsTest, levyConstructor_run )
{
    int N=1;
    int nBeads=10;
    Real beta=1;
    int t0=3;
    int l=4;
    int nTrials=1000;


    SetUp( {N}, nBeads,beta , {300000,300000,300000} );
    SetUpNonInteractingHarmonicAction();


    configurations.setHeadTail(0,t0+l,t0-1);

    pimc::levyReconstructor reconstructor( N * ( nBeads+1) );

    reconstructor.setReconstructorBoundaries( pimc::chainBoundary::free,pimc::chainBoundary::free );
    reconstructor.setSigma(1);
    reconstructor.setGaussianParticleSampling();



    int nBlocks=10000;
    int nStepsPerBlock=10000;


    auto lengthEstimator=std::make_shared< pimc::lengthEstimator >( l );
    lengthEstimator->setStartFromHeadOrTail(0);


    auto lengthObOpen = std::make_shared<pimc::scalarObservable>(lengthEstimator,std::string("lengthOpen"),false);

    configurations.dataTensor().setConstant(0);
    
    std::ofstream x02Out;
    x02Out.open("x02Out.dat");

    const auto & data= configurations.dataTensor();

    for(int iBlock=0;iBlock< nStepsPerBlock;iBlock++)
    {
        Real x02=0;
        for(int n=0;n<nTrials;n++)
        {
            reconstructor.apply( configurations,{t0,t0+l},0,S,randG);
            lengthObOpen->accumulate(configurations,S);
            for(int d=0;d<getDimensions();d++)
            {
                x02+=data(0,d,t0)*data(0,d,t0) ;
            }
        }
    
        lengthObOpen->out(iBlock);
        lengthObOpen->clear();
        x02Out << iBlock << "\t" << x02/nTrials << std::endl << std::flush;
        std::cout << "Block: " << iBlock << std::endl << std::flush;
    }
    x02Out.close();

    // x02: D
    // lengthOpen : l * tau

}


TEST_F( configurationsTest, createRemoveWormSemiCanonical_run )
{

    auto j= R"(
        {
        "CA": 1e-2,
        "CB" : 1e-2,
        "l" : 1,
        "t0" : 0,
        "nBlocks" : 10000,
        "stepsPerBlock" : 100000,
        "randomInitialBead" : false,
        "fixedLength" : false
        }
    )"_json;


    
    if (arguments::testInputFile != "" )
    {
        std::ifstream ifs(arguments::testInputFile);
        j = pimc::json_t::parse(ifs);
    }

    Real CA=j["CA"].get<Real>();
    Real CB=j["CB"].get<Real>();


    int l=j["l"].get<int>();
    int t0=j["t0"].get<int>();
     int nBlocks=j["nBlocks"].get<int>();
    int nStepsPerBlock=j["stepsPerBlock"].get<int>();
    bool randomInitialBead=j["randomInitialBead"].get<bool>();
    bool fixedLength=j["fixedLength"].get<bool>();


    auto jWormMoveSettings = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4,
            "CA" : 1,
            "CB" : 1
        }
            )"_json;

    
    
    jWormMoveSettings["reconstructionMaxLength"]=l;
    jWormMoveSettings["CA"]=CA;
    jWormMoveSettings["CB"]=CB;

    


    
    int nBeads=10;
    int NA=2;
    int NB=0;


    Real beta=0.1* nBeads;
    std::array<Real,3> lBox{ TRUNCATE_D(1,1,1)};
    SetUp( {NA,NB}, nBeads,beta , {300000,300000,300000} );
    SetUpNonInteractingHarmonicAction();
    //SetUpFreeParticleAction();
    SetGrandCanonicalEnsamble( { 0,0});



    const auto & groupA=configurations.getGroups()[0];
    const auto & groupB=configurations.getGroups()[1];


   // SetSeed(456);

    SetRandom(lBox);

    configurations.fillHeads();
    
    int lLong=int(M*0.6);
    int lShort=int(M*0.3);

    int iChain=0;    

    int nBurnsIn=10000;
    
    pimc::createWormSemiCanonicalMove openMove(jWormMoveSettings);
    pimc::removeWormSemiCanonicalMove closeMove(jWormMoveSettings);

    if ( not randomInitialBead)
    {
        openMove.setStartingBead(t0);
        closeMove.setStartingBead(t0);
    }

    if ( fixedLength)
    {
        openMove.setLengthCut(l);
        closeMove.setLengthCut(l);
    }

    pimc::levyMove levyMoveA(lLong,0);
    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);

    pimc::levyMove levyMoveB(lLong,1);
    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);

    Real delta=0.7;
    pimc::translateMove translMoveA(delta,(M+1)*configurations.nChains()*2,0);
    pimc::translateMove translMoveB(delta,(M+1)*configurations.nChains()*2,1);    

    tab.push_back(&levyMoveA,0.8,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::offDiagonal,"translateA");

    tab.push_back(&levyMoveB,0.8,pimc::sector_t::offDiagonal,"levyB");
    tab.push_back(&moveHeadMoveB,0.05,pimc::sector_t::offDiagonal,"moveHeadB");
    tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::offDiagonal,"translateB");


    tab.push_back(&levyMoveA,0.9,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::diagonal,"translateA");

    tab.push_back(&levyMoveB,0.9,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::diagonal,"translateB");

    /* configurations.setHeadTail(groupA.iStart + 0,t0,-1);
    configurations.setHeadTail(groupA.iStart + 1,M,t0 + l -1);
    configurations.join(groupA.iStart + 1,groupA.iStart + 0);

    configurations.setHeadTail(groupB.iStart + 0,t0 + l, t0 - 1);
 */
    
    configurations.fillHeads();

   
   auto lengthEstimatorBLong=std::make_shared< pimc::lengthEstimator >( l );
   auto lengthEstimatorAShort=std::make_shared< pimc::lengthEstimator >( M );
   auto lengthEstimatorALong=std::make_shared< pimc::lengthEstimator >( M - l );

    lengthEstimatorAShort->setStartingBead(groupA.iStart + 0,t0);
    lengthEstimatorBLong->setStartFromHeadOrTail(1);
    lengthEstimatorALong->setStartFromHeadOrTail(0);


    auto lengthAObShort = std::make_shared<pimc::scalarObservable>(lengthEstimatorAShort,std::string("lengthAShort"),false);
    auto lengthAObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorALong,std::string("lengthALong"),false);

    auto lengthBObLong = std::make_shared<pimc::scalarObservable>(lengthEstimatorBLong,std::string("lengthBLong"),false);

    std::cout << "Equilibrating..." << std::endl;

    for(int n=0;n<nBurnsIn;n++)
    {
        tab.attemptMove(configurations,S,randG);
    }

    const auto & data = configurations.dataTensor();

    Real nShort=0;
    Real nLong=0;
    Real ratio=0;
    Real n=0;
    

    std::ofstream openRatioOut("openRatio.dat");

    std::cout << "Start sampling..." << std::endl;
    for(int t=0;t<nBlocks;t++)
    {
       while ( (nShort < nStepsPerBlock) and (nLong < nStepsPerBlock) and ( n < nStepsPerBlock ) )
        {
            n+=1;
            bool accept=tab.attemptMove(configurations,S,randG);

            if ( not configurations.isOpen(0) )
            {
                openMove.attemptMove(configurations,S,randG);
            }
            else 
            {
                closeMove.attemptMove(configurations,S,randG);
            }

            if ( not configurations.isOpen(0) )
            {
                nShort++;
                ratio+=1;
                lengthAObShort->accumulate(configurations,S);
            }
            else
            {
                lengthAObLong->accumulate(configurations,S);
                lengthBObLong->accumulate(configurations,S);
                nLong++;
            }
        }

    
    std::cout << "Short length fraction: " << nShort*1./n << std::endl;

    tab >> std::cout ;

    if (nShort == nStepsPerBlock)
        {
           

            lengthAObShort->out(t);
            lengthAObShort->clear();

            nShort=0;
        }
    
    if ( nLong== nStepsPerBlock)
        {
            lengthBObLong->out(t);
            lengthBObLong->clear();

            lengthAObLong->out(t);
            lengthAObLong->clear();

            nLong=0;
        }

    if ( n == nStepsPerBlock)
        {
            ratio/=1.*n;
            openRatioOut << t << " " << ratio << std::endl;
            ratio=0;
            n=0;
        }

    }
}



bool checkHeadTailMatch( const pimc::configurations_t & confs, int setA, int setB)
    {
        if (confs.isOpen(setA) and confs.isOpen(setB) ) 
        {
            auto iChainHeadA = confs.getGroups()[setA].heads[0];
            auto iChainTailB = confs.getGroups()[setB].tails[0];

            auto iChainHeadB = confs.getGroups()[setB].heads[0];
            auto iChainTailA = confs.getGroups()[setA].tails[0];

            auto tHeadA=confs.getChain(iChainHeadA).head;
            auto tTailB=confs.getChain(iChainTailB).tail;

            auto tHeadB=confs.getChain(iChainHeadB).head;
            auto tTailA=confs.getChain(iChainTailA).tail;

            return (tHeadA == (tTailB + 1)) & (tHeadB == (tTailA + 1));

        }

        return true;
    }

#include "../pimc/exampleActions.h"

TEST_F( configurationsTest, semiCanonical_run )
{

    auto j= R"(
        {
        "CA": 1e-2,
        "CB" : 1e-2,
        "lLong" : 6,
        "lShort" : 4,
        "nBlocks" : 10000,
        "stepsPerBlock" : 100000,
        "NA" : 2,
        "NB" : 0,
        "deltaMu" : 0
        }
    )"_json;


    SetSeed(236789);




    if (arguments::testInputFile != "" )
    {
        std::ifstream ifs(arguments::testInputFile);
        j = pimc::json_t::parse(ifs);
    }

    int NA=j["NA"].get<int>();
    int NB=j["NB"].get<int>();

    Real CA=j["CA"].get<Real>();
    Real CB=j["CB"].get<Real>();

    Real deltaMu= j["deltaMu"].get<Real>();

    int lLong=j["lLong"].get<int>();
    int lShort=j["lShort"].get<int>();
    
    int nBlocks=j["nBlocks"].get<int>();
    int nStepsPerBlock=j["stepsPerBlock"].get<int>();
    
    auto jWormMoveSettingsA = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4,
            "CA" : 1,
            "CB" : 1,
            "initialBeadSampling" : {
                "kind" : "uniform"
            }
        }
            )"_json;

     auto jAdvanceHeadTailA = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4
        }
            )"_json;
     
     auto jAdvanceHeadTailB = R"(
        {
            "setA": 1,
            "setB": 0 ,
            "reconstructionMaxLength": 4
        }
            )"_json;

    
    jAdvanceHeadTailB["reconstructionMaxLength"]=lShort;
    jAdvanceHeadTailA["reconstructionMaxLength"]=lShort;
            
    
    jWormMoveSettingsA["reconstructionMaxLength"]=lShort;
    jWormMoveSettingsA["CA"]=CA;
    jWormMoveSettingsA["CB"]=CB;

    auto jWormMoveSettingsB = jWormMoveSettingsA;
    jWormMoveSettingsB["CA"]=CB;
    jWormMoveSettingsB["CB"]=CA;
    jWormMoveSettingsB["setA"]=1;
    jWormMoveSettingsB["setB"]=0;

     auto jSemiOpenA = R"(
        {
            "setA": 0,
            "setB": 1,
            "reconstructionMaxLength": 4,
            "C" : 1
        }
            )"_json;


     auto jSemiOpenB = R"(
        {
            "setA": 1,
            "setB": 0,
            "reconstructionMaxLength": 2,
            "C" : 1
        }
            )"_json;
    
    

    jSemiOpenA["C"]=CA;
    jSemiOpenA["reconstructionMaxLength"]=lShort;


    jSemiOpenB["C"]=CB;
    jSemiOpenB["reconstructionMaxLength"]=lShort;



    pimc::semiOpenMove semiOpenA(jSemiOpenA);
    pimc::semiCloseMove semiCloseA(jSemiOpenA);

    pimc::semiOpenMove semiOpenB(jSemiOpenB);
    pimc::semiCloseMove semiCloseB(jSemiOpenB);

    pimc::fullSemiCanonicalOpenMove openFullMoveA(jSemiOpenA);
    pimc::fullSemiCanonicalCloseMove closeFullMoveA(jSemiOpenA);

    pimc::fullSemiCanonicalOpenMove openFullMoveB(jSemiOpenB);
    pimc::fullSemiCanonicalCloseMove closeFullMoveB(jSemiOpenB);

    pimc::swapMove swapMoveA( lShort, NA + NB , 0 );
    pimc::swapMove swapMoveB( lShort, NA + NB , 1 );

    //swapMoveA.setFixedLength();

    //openFullMoveA.setLength(lShort);
    //closeFullMoveA.setLength(lShort);


    //openFullMoveB.setLength(lShort);
    //closeFullMoveB.setLength(lShort);

    int nBeads=10;

    Real beta=1;


    std::array<Real,3> lBox{ TRUNCATE_D(10,10,10)};
    SetUp( {NA,NB}, nBeads,beta , lBox );



    //SetUpNonInteractingHarmonicAction();
    SetUpFreeParticleAction();
    //SetUpTwoBodyInteractionHarmonicInTrap( { {0,1 } , {0,0}, {1,1} });

    Real a=0.1;
    auto S = createFreeCaoBerneAction( a , timeStep ,  geo,{0,0}, configurations.nChains() , configurations.nBeads() );
    //SetUpFreeParticleAction();
    SetGrandCanonicalEnsamble( {deltaMu,0});

    const auto & groupA=configurations.getGroups()[0];
    const auto & groupB=configurations.getGroups()[1];

    
    SetRandom(lBox);

    configurations.fillHeads();

    int iChain=0;
    int nBurnsIn=10000;

    //semiOpenA.setStartingBead(5);
    //semiCloseA.setStartingBead(5);

    //semiOpenA.setStartingChain(groupA.iStart);
    //semiCloseA.setStartingChain(groupA.iStart);
    //semiCloseA.setLength(lShort);



    //semiOpenB.setStartingBead(5);
    //semiCloseB.setStartingBead(5);
    //semiOpenB.setStartingChain(groupB.iStart);
    //semiCloseB.setStartingChain(groupB.iStart);


    //openFullMoveB.setStartingChain(groupB.iStart);
    //closeFullMoveB.setStartingChain(groupB.iStart);

    //openFullMoveA.setStartingChain(groupA.iStart);
    //closeFullMoveA.setStartingChain(groupA.iStart);


    

    pimc::createWormSemiCanonicalMove openMoveA(jWormMoveSettingsA);
    pimc::removeWormSemiCanonicalMove closeMoveA(jWormMoveSettingsA);

    pimc::createWormSemiCanonicalMove openMoveB(jWormMoveSettingsB);
    pimc::removeWormSemiCanonicalMove closeMoveB(jWormMoveSettingsB);




    //openMoveA.setInitialUniformSampling();
    //openMoveA.setInitialGaussianSampling(1);

    //closeMoveA.setInitialUniformSampling();
    //closeMoveA.setInitialGaussianSampling(1);

    //openMoveB.setInitialUniformSampling();
    //openMoveB.setInitialGaussianSampling(1);

    //closeMoveB.setInitialUniformSampling();
    //openMoveB.setInitialGaussianSampling(1);

    //openMoveA.setLengthCut( lShort);
    //closeMoveA.setLengthCut( lShort);

    //openMoveA.setStartingBead( 6);
    //closeMoveA.setStartingBead( 6);
    

  



    pimc::advanceHeadTail advanceMoveA(jAdvanceHeadTailA);
    pimc::recedeHeadTail recedeMoveA(jAdvanceHeadTailA);

    pimc::advanceHeadTail advanceMoveB(jAdvanceHeadTailB);
    pimc::recedeHeadTail recedeMoveB(jAdvanceHeadTailB);


    pimc::levyMove levyMoveA(lLong,0);
    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);

    pimc::levyMove levyMoveB(lLong,1);
    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);

    Real delta=10;
    pimc::translateMove translMoveA(delta,(M+1)*configurations.nChains()*2,0);
    pimc::translateMove translMoveB(delta,(M+1)*configurations.nChains()*2,1);    


    std::shared_ptr<pimc::scalarEstimator> nAEst=std::make_shared<pimc::particleNumberEstimator>( 0 );
    std::shared_ptr<pimc::scalarEstimator>  nBEst=std::make_shared<pimc::particleNumberEstimator>( 1 );
    std::shared_ptr<pimc::scalarEstimator>  nEst=std::make_shared<pimc::particleNumberEstimator>( std::vector<int>{0,1} );


    auto pairEst= std::make_shared<pimc::pairCorrelation>( 0,0 );



    auto nRingsAEst = std::make_shared<pimc::nConnectedChainsEstimator>(0);
    auto nRingsBEst = std::make_shared<pimc::nConnectedChainsEstimator>(1);

    auto totalLengthEstimatorA = std::make_shared<pimc::closedLengthEstimator>(0);




    auto nAOb = std::make_shared< pimc::scalarObservable>(nAEst,"NA",false);
    auto nBOb = std::make_shared< pimc::scalarObservable>(nBEst,"NB",false);
    auto nOb = std::make_shared< pimc::scalarObservable>(nEst,"N",false);
    auto pairOb = std::make_shared< pimc::histogramObservable>(pairEst,"pair",100,0,5);

    
    auto nRingsAOb = std::make_shared< pimc::scalarObservable>(nRingsAEst,"nRingsA",false);
    auto nRingsBOb = std::make_shared< pimc::scalarObservable>(nRingsBEst,"nRingsB",false);
     auto totalLengthAOb = std::make_shared< pimc::scalarObservable>(totalLengthEstimatorA,"totalLengthA",false);



    auto wormLengthAEst=std::make_shared<pimc::nBeadsInWormEstimator>( 0 );
    auto wormLengthBEst=std::make_shared<pimc::nBeadsInWormEstimator>( 1 );

    auto wormLengthAOb = std::make_shared< pimc::scalarObservable>(wormLengthAEst,"lWormA",false);
    auto wormLengthBOb = std::make_shared< pimc::scalarObservable>(wormLengthBEst,"lWormB",false);

    auto lengthEstimatorA=std::make_shared< pimc::lengthEstimator >(  );
    auto lengthEstimatorB=std::make_shared< pimc::lengthEstimator >( );

    auto l2WormA=std::make_shared< pimc::lengthEstimator >(  );
    l2WormA->setStartFromHeadOrTail(0);


    lengthEstimatorA->setStartingBead(groupA.iStart + 0,0);
    lengthEstimatorB->setStartingBead(groupB.iStart + 0,0);


    auto l2WormAOb = std::make_shared<pimc::scalarObservable>(l2WormA,std::string("l2WormA"),false);
    

    auto lengthAOb = std::make_shared<pimc::scalarObservable>(lengthEstimatorA,std::string("lengthChainA"),false);
    auto lengthBOb = std::make_shared<pimc::scalarObservable>(lengthEstimatorB,std::string("lengthChainB"),false);


    

    auto eEstim=std::make_shared<pimc::thermodynamicEnergyEstimator> ();
    auto eVEstim = std::make_shared<pimc::virialEnergyEstimator> (configurations.nChains(),configurations.nBeads() );

    auto eOb = std::make_shared<pimc::scalarObservable>(eEstim,std::string("e"),false);
    auto eVOb = std::make_shared<pimc::scalarObservable>(eVEstim,std::string("eV"),false);


    tab.push_back(&levyMoveA,0.75,pimc::sector_t::offDiagonal,"levyA");
    tab.push_back(&moveHeadMoveA,0.05,pimc::sector_t::offDiagonal,"moveHeadA");
    tab.push_back(&moveTailMoveA,0.05,pimc::sector_t::offDiagonal,"moveTailA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::offDiagonal,"translateA");
    //tab.push_back(&advanceMoveA,0.05,pimc::sector_t::offDiagonal,"advanceMoveA");
    //tab.push_back(&recedeMoveA,0.05,pimc::sector_t::offDiagonal,"recedeMoveA");
    //tab.push_back(&closeMoveA,0.05,pimc::sector_t::offDiagonal,"closeMoveA");
    tab.push_back(&semiCloseA,0.05,pimc::sector_t::offDiagonal,"semiCloseMoveA");
    //tab.push_back(&closeFullMoveA,0.05,pimc::sector_t::offDiagonal,"closeFullMoveA");
    //tab.push_back(&swapMoveA,0.1,pimc::sector_t::offDiagonal,"swapMoveA");


    //tab.push_back(&levyMoveB,0.5,pimc::sector_t::offDiagonal,"levyB");
    //tab.push_back(&moveHeadMoveB,0.05,pimc::sector_t::offDiagonal,"moveHeadB");
    //tab.push_back(&moveTailMoveB,0.05,pimc::sector_t::offDiagonal,"moveTailB");
    //tab.push_back(&translMoveB,0.05,pimc::sector_t::offDiagonal,"translateB");
    //tab.push_back(&advanceMoveB,0.05,pimc::sector_t::offDiagonal,"advanceMoveB");
    //tab.push_back(&recedeMoveB,0.05,pimc::sector_t::offDiagonal,"recedeMoveB");

    //tab.push_back(&closeMoveB,0.05,pimc::sector_t::offDiagonal,"closeMoveB");
    //tab.push_back(&closeFullMoveB,0.05,pimc::sector_t::offDiagonal,"closeFullMoveB");
    //tab.push_back(&semiCloseB,0.05,pimc::sector_t::offDiagonal,"semiCloseMoveB");
    //tab.push_back(&swapMoveB,0.1,pimc::sector_t::offDiagonal,"swapMoveB");

    tab.push_back(&levyMoveA,0.85,pimc::sector_t::diagonal,"levyA");
    tab.push_back(&translMoveA,0.1,pimc::sector_t::diagonal,"translateA");
    //tab.push_back(&openMoveA,0.05,pimc::sector_t::diagonal,"openMoveA");
    tab.push_back(&semiOpenA,0.05,pimc::sector_t::diagonal,"semiOpenMoveA");
    //tab.push_back(&openFullMoveA,0.05,pimc::sector_t::diagonal,"openFullMoveA");

    tab.push_back(&levyMoveB,0.9,pimc::sector_t::diagonal,"levyB");
    tab.push_back(&translMoveB,0.1,pimc::sector_t::diagonal,"translateB");
    //tab.push_back(&openMoveB,0.05,pimc::sector_t::diagonal,"openMoveB");
    //tab.push_back(&openFullMoveB,0.05,pimc::sector_t::diagonal,"openFullMoveB");
    //tab.push_back(&semiOpenB,0.05,pimc::sector_t::diagonal,"semiOpenMoveB");


    /*  int t0=3;
    configurations.setHeadTail(0 + groupA.iStart,  M , t0-1);
    configurations.setHeadTail(1 + groupA.iStart,  t0  , -1);
    configurations.join(0 + groupA.iStart,1 + groupA.iStart);
 */
  /*   configurations.setHeadTail(0 + groupB.iStart,  M , t0  -1);
    configurations.setHeadTail(1 + groupB.iStart,  t0 , -1);
    configurations.join(0 + groupB.iStart,1 + groupB.iStart); */

  /*   int t0=5;
    configurations.setTail(groupA.iStart + 0,t0-1);
    configurations.setHead(groupA.iStart + 1,t0);
    configurations.join(groupA.iStart + 0,groupA.iStart + 1);
     */
    
    //configurations.join(groupA.iStart + 0,groupA.iStart + 1);
    //configurations.join(groupA.iStart + 1,groupA.iStart + 0);

    //configurations.setTail(groupB.iStart + 0,t0-1);
    //configurations.setHead(groupB.iStart + 1,t0);
    //configurations.join(groupB.iStart + 0,groupB.iStart + 1);


    //int t0=3;
    //configurations.setHead(groupA.iStart + 2, t0);
    //configurations.setTail(groupA.iStart + 1, t0-1);
    //configurations.join(groupA.iStart + 1 ,groupA.iStart + 2);

    //configurations.join(groupA.iStart + 2 ,groupA.iStart + 0);


    configurations.fillHeads();
    {
        int iAttemptInitialCondition=0;

        while (not S.checkConstraints(configurations) )
        {
            configurations.setRandom( { geo.getLBox(0) ,geo.getLBox(1),geo.getLBox(2) } , randG );
            iAttemptInitialCondition++;

            if (iAttemptInitialCondition > 100000)
            {
                throw std::runtime_error("Max iteration reached in generating the initial condition");
            }
        }
    }
/* 
    std::cout << "Equilibrating..." << std::endl;

    for(int n=0;n<nBurnsIn;n++)
    {
        tab.attemptMove(configurations,S,randG);
    } */

    const auto & data = configurations.dataTensor();

    Real nShort=0;
    Real nLongA=0;
    Real nLongB=0;

    Real ratio=0,ratioA=0,ratioB=0;
    Real n=0;


    std::ofstream closedRatioOut("closedRatio.dat");
    std::ofstream closedRatioAOut("closedRatioA.dat");
    std::ofstream closedRatioBOut("closedRatioB.dat");
    
    

     if ( ! fs::exists("configurations") ) 
    { 
        fs::create_directory("configurations"); // create src folder
    }

    std::cout << "Start sampling..." << std::endl;
    for(int t=0;t<nBlocks;t++)
    {
       while ( (nShort < nStepsPerBlock) and (nLongA < nStepsPerBlock) and ( n < nStepsPerBlock )  and (nLongB < nStepsPerBlock) )
        {
            n+=1;
            bool accept=tab.attemptMove(configurations,S,randG);

            //bool valid = checkTimePeriodicBoundaryConditions(configurations,geo) and checkHeadTailMatch( configurations, 0, 1);

            if ( not configurations.isOpen() )
            {
                //auto N = configurations.nParticles();


                nShort++;
                nAOb->accumulate(configurations,S);
                nBOb->accumulate(configurations,S);
                nOb->accumulate(configurations,S);
                nRingsAOb->accumulate(configurations,S);
                nRingsBOb->accumulate(configurations,S);

                lengthAOb->accumulate(configurations,S);
                lengthBOb->accumulate(configurations,S);

                eOb->accumulate(configurations,S);
                eVOb->accumulate(configurations,S);

                totalLengthAOb->accumulate(configurations,S);
                pairOb->accumulate(configurations,S);

                ratio+=1;
                ratioA+=1;
                ratioB+=1;

            }
            else
            {   

                if (configurations.isOpen(0)   )
                {
                    nLongA++;
                    wormLengthAOb->accumulate(configurations,S);
                    l2WormAOb->accumulate(configurations,S);

                    
                } else
                {
                    ratioA+=1;

                }
                if (configurations.isOpen(1))
                {
                    nLongB++;
                    wormLengthBOb->accumulate(configurations,S);
                }
                else
                {
                    ratioB+=1;
                }
                
            }
        }

    std::cout << "Short length fraction: "<< nShort*1./n << std::endl;

    tab >> std::cout ;

    if (nShort == nStepsPerBlock)
        {
           nAOb->out(t);
           nAOb->clear();

           nBOb->out(t);
           nBOb->clear();

           nOb->out(t);
           nOb->clear();

           nRingsBOb->out(t);
           nRingsBOb->clear();

           nRingsAOb->out(t);
           nRingsAOb->clear();

           lengthAOb->out(t);
           lengthAOb->clear();

           lengthBOb->out(t);
           lengthBOb->clear();

           eOb->out(t);
           eOb->clear();

           eVOb->out(t);
           eVOb->clear();

           totalLengthAOb->out(t);
           totalLengthAOb->clear();

           pairOb->out(t);
           pairOb->clear();

           

           nShort=0;
        }
    
    if ( nLongA== nStepsPerBlock)
        {
            wormLengthAOb->out(t);
            wormLengthAOb->clear();

            l2WormAOb->out(t);
            l2WormAOb->clear();

            nLongA=0;
        }
    if ( nLongB== nStepsPerBlock)
        {
            wormLengthBOb->out(t);
            wormLengthBOb->clear();

            
            

            nLongB=0;
        }


    if ( n == nStepsPerBlock)
        {
            ratio/=1.*n;
            ratioA/=1.*n;
            ratioB/=1.*n;

            closedRatioOut << t << " " << ratio << std::endl;
            closedRatioAOut << t << " " << ratioA << std::endl;
            closedRatioBOut << t << " " << ratioB << std::endl;

            ratio=0;ratioA=0;ratioB=0;
            n=0;
        }

             configurations.saveHDF5("configurations/particles" + std::to_string(t) + ".hdf5");
        
    }

  

}


TEST_F( configurationsTest, restrictInBox )
{
    std::array<Real,getDimensions() > lBox { 1, 1, 1};
    Real beta=1;

    int NA=20;
    int NB=30;
    int nBeads=100;

    SetUp( {NA,NB}, nBeads,beta , lBox );

    configurations.setRandom( { lBox[0] *5 ,lBox[1]*5,lBox[2]*5 } , randG );
    Eigen::Tensor<Real,3> dataOld(configurations.dataTensor());
    
    restrictToBox(configurations,geo);

    const auto & data = configurations.dataTensor();

    const auto & groups=configurations.getGroups();


    for (const auto & group : groups)
    {
        for(int i=group.iStart;i<=group.iEnd;i++)
        {
            const auto & chain = configurations.getChain(i);

            for(int d=0;d<getDimensions();d++)
            {    
                ASSERT_LE( data(i,d,chain.tail+1) , lBox[d]/2);
                ASSERT_GE( data(i,d,chain.tail+1) , -lBox[d]/2);
            }

            for(int t=0;t<M;t++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    auto newD=data(i,d,t+1) - data(i,d,t);
                    auto oldD=dataOld(i,d,t+1) - dataOld(i,d,t);
                    ASSERT_NEAR(newD,oldD, TOL);

                }
                
            }
        }
    }

}