#include "integrationTests.h"

#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"


void harmonicTrapTest::SetUpTwoBodyInteractionHarmonicInTrap_kernel()
{
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    int order=2;

    auto V = pimc::makeIsotropicPotentialFunctor(
        [](Real r) {return 0.5*(r*r) ;} ,
        [](Real r) {return r  ;} );

        auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);


    auto sV2B=std::make_shared<pimc::actionTwoBody>();
    sV2B->setSets({0,0});
    sV2B->setKernel(kernel);


    auto  sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,order);

    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,sV2B};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

}

void twoBodyTest::SetUpFreeActionWithHardSphereConstraint(Real a)
{
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    auto V = pimc::makeIsotropicPotentialFunctor(
        [](Real r) {return 0 ;} ,
        [](Real r) {return 0  ;} );

        auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);


    auto sV2B=std::make_shared<pimc::actionTwoBody>();
    sV2B->setSets({0,0});
    sV2B->setKernel(kernel);
    sV2B->setGeometry(geo);
    sV2B->setTimeStep(timeStep);
    sV2B->setMinimumDistance(a);
    
    std::shared_ptr<pimc::action> sEmpty= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);

    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sV2B};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

}

void twoBodyTest::SetUpCaoBernePropagator(Real radius)
{
    using propagator_t = pimc::caoBernePropagator;

    auto kernel =std::make_shared<propagator_t>(timeStep,radius);
    auto pairKernel = std::make_shared<pimc::pairProductKernel<propagator_t> >(kernel);

    auto S2B = std::make_shared<pimc::actionTwoBody>();
    S2B->setTimeStep(timeStep);
    S2B->setGeometry(geo);
    S2B->setKernel(pairKernel);
    pairKernel->setGeometry(geo);
    pairKernel->setTimeStep(timeStep);

    S2B->setSets({0,0});
    S2B->setMinimumDistance(radius);
    S2B->setGeometry(geo);
    S2B->setTimeStep(timeStep);


    auto  sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , configurations.nChains()  , geo);

    std::shared_ptr<pimc::action> sEmpty= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);

    


    std::vector<std::shared_ptr<pimc::action> > Vs = {S2B};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

}


void twoBodyTest::SetRandomMinimumDistance(Real radius, const std::array<Real,DIMENSIONS> lBox)
{

    while (not S.checkConstraints(configurations) )
    {
        SetRandom(lBox);
    } ;
}


void twoBodyTest::sample()
{

    Real C=1e-3;
    int nBeads=10;
    int N=1;
    Real beta=0.1* nBeads;

    SetUp(N,nBeads,beta );

    //SetUpFreeParticleAction();    
    //SetUpNonInteractingHarmonicAction();

    //SetUpTwoBodyInteractionHarmonic();
    SetUpTwoBodyInteractionHarmonicInTrap();

    SetGrandCanonicalEnsamble( 2);
    SetSeed( time(NULL) );
    SetRandom({TRUNCATE_D(0.4,0.4,0.4)});

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
    int lOpen=3;

    pimc::translateMove translate(0.1, 2000*M , 0 );

    pimc::levyMove levy(l,0);
    
    pimc::moveHead moveHeadMove(lShort,0);
    pimc::moveTail moveTailMove(lShort,0);

    pimc::openMove open(C, 0, lOpen );
    pimc::closeMove close(C, 0, lOpen );
    
    pimc::createWorm createWorm(C, 0, lShort , 1 );
    pimc::deleteWorm removeWorm(C, 0, lShort , 1);
    

    //open.setStartingBead(7);
    //open.setStartingChain(0);

    //close.setStartingBead(7);
    //close.setStartingChain(0);

    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);
    
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


void twoBodyTest::accumulate()
{
    resetCounters();

    
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
    pimc::nConnectedChains nConnectedChains;

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
                l2Long+=accumulateAverageLengthSquare( 0,configurations );
                //l2Long+=accumulateLengthSquare( configurations , {0,N-1} , {0,nBeads-1} , geo );

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


TEST_F(twoBodyTest,pbc)
{
    Real C=1e-2;
    int nBeads=10;
    int N=2;
    Real beta=0.1* nBeads;

    SetUp(N,nBeads,beta , {1} );
    nTrials*=10;

    SetUpFreeParticleAction();
    SetGrandCanonicalEnsamble( -1 );
    SetSeed( 145) ;
    SetRandom({TRUNCATE_D(0.4,0.4,0.4)});

    int t0=7;
    int l = int( 0.8* 10);
    int lShort=int( 0.6* 10);
    int lOpen=1;

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

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();

    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);
    //swap.setFixedLength();
    
    //advanceHead.setMaximumParticleNumber(2);
    //recedeHead.setMinParticleNumber(1);

    pimc::nConnectedChains nConnectedChains;
    tab.push_back(&levy,0.5,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.4,pimc::sector_t::diagonal,"translate");
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
    
    int tTail=4;
    int lWormShort=4;

    //configurations.setHeadTail(0,tTail + lWormShort,tTail -1);

    //configurations.setHeadTail(1,tTail - lCut,-1);
    //configurations.join(0,1);



    //configurations.join(1,0);
        

    configurations.fillHeads();

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;

    accumulate();

}


TEST_F(harmonicTrapTest, twoBodyActionKernel2_grandCanonical )
{
    Real C=1;
    int nBeads=10;
    int N=2;
    Real beta=0.1* nBeads;

    SetUp(N,nBeads,beta , {10000000} );
    //SetUpTwoBodyInteractionHarmonicInTrap();
    SetUpTwoBodyInteractionHarmonicInTrap_kernel();

    SetGrandCanonicalEnsamble( 0 );
    SetSeed( 145) ;
    SetRandom({TRUNCATE_D(0.4,0.4,0.4)});

    int t0=7;
    int l = int( 0.8* 10);
    int lShort=int( 0.6* 10);
    int lOpen=1;

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

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();

    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);
    //swap.setFixedLength();
    
    //advanceHead.setMaximumParticleNumber(2);
    //recedeHead.setMinParticleNumber(1);


    pimc::nConnectedChains nConnectedChains;
    tab.push_back(&levy,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.4,pimc::sector_t::diagonal,"translate");
    //tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");

    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");

    tab.push_back(&levy,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&advanceHead,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.05,pimc::sector_t::offDiagonal, "recedeHead");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");

    int tTail=4;
    int lWormShort=4;

    //configurations.setHeadTail(0,tTail + lWormShort,tTail -1);

    //configurations.setHeadTail(1,tTail - lCut,-1);
    //configurations.join(0,1);

    //configurations.join(1,0);        

    configurations.fillHeads();

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;

    accumulate();
    
}


#if DIMENSIONS == 3
TEST_F( twoBodyTest, caoBernePropagator )
{
    Real C=1;
    int nBeads=10;

    int N=2;
    Real beta=0.1* nBeads;
    Real a=0.1;

    std::array<double,DIMENSIONS> lBox = {TRUNCATE_D(1,1,1)};

    SetUp(N,nBeads,beta , lBox );
    //SetUpTwoBodyInteractionHarmonicInTrap();
    SetUpCaoBernePropagator(a);
    //SetUpFreeParticleAction();


    //SetUpFreeActionWithHardSphereConstraint(a);    
    SetGrandCanonicalEnsamble( 0 );
    SetSeed( time( NULL)  ) ;
    SetRandomMinimumDistance(a,lBox);


    int t0=7;
    int l = int( 0.3* 10);
    int lShort=int( 0.6* 10);
    int lOpen=1;

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

    pimc::advanceHead advanceHead(lShort,0);
    pimc::recedeHead recedeHead(lShort,0);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();

    pimc::advanceTail advanceTail(lShort,0);
    pimc::recedeTail recedeTail(lShort,0);

    pimc::swapMove swap( lShort , 200 , 0);
    //swap.setFixedLength();
    
    //advanceHead.setMaximumParticleNumber(2);
    //recedeHead.setMinParticleNumber(1);

    pimc::nConnectedChains nConnectedChains;
    tab.push_back(&levy,0.9,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::diagonal,"translate");
    //tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");

    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");


    tab.push_back(&levy,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&advanceHead,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHead,0.05,pimc::sector_t::offDiagonal, "recedeHead");
    tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");
    
    int tTail=4;
    int lWormShort=4;

    //configurations.setHeadTail(0,tTail + lWormShort,tTail -1);

    //configurations.setHeadTail(1,tTail - lCut,-1);
    //configurations.join(0,1);
    //configurations.join(1,0);        

    configurations.fillHeads();
    resetCounters();

    int nTrials=10000;
    int nBlocks=100000;

    auto g = std::make_shared<pimc::pairCorrelation>(0,0);
    auto gOb=std::make_shared<pimc::histogramObservable>(g,"pairCorr",100,0,lBox[0]/2 );
    
    auto eVEst = std::make_shared<pimc::virialEnergyEstimator>(configurations.nChains() , configurations.nBeads()  );

    auto eEst= std::make_shared< pimc::thermodynamicEnergyEstimator>();
    
    auto eO= std::make_shared<pimc::scalarObservable>(eEst,std::string("energy") );

    auto eVO= std::make_shared<pimc::scalarObservable>(eVEst,std::string("eV") );


    int nUnCorrelationSteps = 100;

    configurations.save( "configurationsInitial" , "csv");

    for (int i=0;i<nBlocks;i++)
    {
        int nClosed=0;
        int nOpen=0;
        int n=0;

        while (nClosed<nTrials and nOpen<nTrials and n<nTrials)
        {
            for(int n=0;n<nUnCorrelationSteps;n++)
            {
                tab.attemptMove(configurations,S,randG);
            }
            
            if ( configurations.isOpen() )
            {
                nOpen+=1;
            }
            else 
            {
                gOb->accumulate(configurations,S);
                nClosed+=1;
                eVO->accumulate(configurations,S);
                eO->accumulate(configurations,S);
            }

            n+=1;
        }


        ASSERT_TRUE(S.checkConstraints(configurations));


       if ( nClosed == nTrials  )
       {
        gOb->out(i);
        gOb->clear();

        eO->out(i);
        eO->clear();

        eVO->out(i);
        eVO->clear();

       }
       

       tab >> std::cout;
       tab.resetCounters();
       configurations.save( "configurations" + std::to_string(i) , "csv");

    }
    
}

#endif