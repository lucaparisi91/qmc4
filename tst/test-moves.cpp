#include "integrationTests.h"

#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"
#include "testConfigurations.h"


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
            "group": 0,
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
            "group": 0,
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
    SetGrandCanonicalEnsamble( 0 );


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

    configurations.setHeadTail(0 + groupA.iStart,  M , t0-1);
    configurations.setHeadTail(1 + groupA.iStart,  t0  , -1);
    configurations.join(0 + groupA.iStart,1 + groupA.iStart);

    configurations.setHeadTail(0 + groupB.iStart,  M , t0  -1);
    configurations.setHeadTail(1 + groupB.iStart,  t0 , -1);
    configurations.join(0 + groupB.iStart,1 + groupB.iStart);

    //configurations.setHeadTail(1,  t0 , -1);
    //configurations.join(0,1); */

    configurations.fillHeads();


   auto lengthEstimatorAShort=std::make_shared< pimc::lengthEstimator >( M );
   auto lengthEstimatorALong=std::make_shared< pimc::lengthEstimator >( M + l );
   auto lengthEstimatorBShort=std::make_shared< pimc::lengthEstimator >( M );
   auto lengthEstimatorBLong=std::make_shared< pimc::lengthEstimator >( M - l );


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

