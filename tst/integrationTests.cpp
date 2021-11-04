#include "integrationTests.h"

#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"

void harmonicTrapTest::SetUpTwoBodyInteractionHarmonicInTrap_kernel()
{
    SetUpTwoBodyInteractionHarmonicInTrap_kernel({ {0,0} } );
}

void harmonicTrapTest::SetUpTwoBodyInteractionHarmonicInTrap_kernel(const std::vector<std::pair<int,int> > & sets )
{
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    int order=2;

    auto V = pimc::makeIsotropicPotentialFunctor(
        [](Real r) {return 0.5*(r*r) ;} ,
        [](Real r) {return r  ;} );

        auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);

    std::vector<std::shared_ptr<pimc::action> > Vs;

    for ( auto [ setA, setB ] : sets )
        {
            auto sV2B=std::make_shared<pimc::actionTwoBody>();
            sV2B->setSets({setA,setB});
            sV2B->setKernel(kernel);
            sV2B->setTimeStep(timeStep);
            sV2B->setGeometry(geo);
            Vs.push_back(sV2B);

        }

    for (int set=0;set<configurations.getGroups().size() ; set++)
    {
        auto  sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,set,order);
        Vs.push_back(sOneBody);

    }


    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

}


void twoBodyTest::SetUpTwoBodyInteractionGaussian_kernel(Real V0, Real R0)
{
    std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);
    
    Real alpha=1./(2*R0*R0);

    auto V = pimc::makeIsotropicPotentialFunctor(
        [V0, alpha](Real r) {return V0*exp(-alpha*r*r) ;} ,
        [V0,alpha](Real r) {return  V0*(-2*r*alpha)*exp(-alpha*r*r) ;} );

    auto kernel = std::make_shared<                         pimc::primitiveApproximationTwoBodyKernel<decltype(V) >  >(std::make_shared<decltype(V)>(V));

    kernel->setTimeStep(timeStep);
    kernel->setGeometry(geo);


    auto sV2B=std::make_shared<pimc::actionTwoBody>();
    sV2B->setSets({0,0});
    sV2B->setKernel(kernel);
    sV2B->setTimeStep(timeStep);
    sV2B->setGeometry(geo);


    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sV2B};

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

#if DIMENSIONS == 3

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


void twoBodyTest::SetUpCaoBernePropagatorTrapped(Real radius,Real omega)
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
    int order=2;

    auto  sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , configurations.nChains()  , geo);

    std::shared_ptr<pimc::action> sEmpty= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);

    auto V = pimc::makeIsotropicPotentialFunctor(
         [omega](Real r) {return omega*(r*r) ;} ,
         [omega](Real r) {return 2*omega*r  ;} );
    
     auto  sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,order);



    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,S2B};


    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

}



#endif

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

    SetGrandCanonicalEnsamble( 2  );
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

    int nTrials=1000;
    int nBlocks=100000;

    accumulate();

}

TEST_F( harmonicTrapTest ,mixture_twoBody)
{
    Real C=1e-3;
    int nBeads=10;
    int N1=2;
    int N2=2;


    Real beta= 0.1* nBeads;

    SetUp( {N1,N2} ,nBeads,beta , {TRUNCATE_D(10000,10000,10000)} );

    //SetUpFreeParticleAction();    
    
    //SetUpNonInteractingHarmonicAction();


    //SetUpTwoBodyInteractionHarmonic();
    
    SetUpTwoBodyInteractionHarmonicInTrap_kernel( { {0,0} ,{0,1} , {1,1}} );


    //SetUpTwoBodyInteractionHarmonicInTrap_kernel( { } );


    SetGrandCanonicalEnsamble( {0,0} );
    SetSeed( time(NULL) );
    SetRandom({TRUNCATE_D(1,1,1)});

    int l = int( 0.8* 10);
    int lShort=int( 0.3* 10);
    int lOpen=3;


    pimc::translateMove translateA(0.1, 2000*M , 0 );
    pimc::translateMove translateB(0.1, 2000*M , 1 );


    pimc::levyMove levyA(l,0);
    pimc::levyMove levyB(l,1);
    

    pimc::moveHead moveHeadMoveA(lShort,0);
    pimc::moveTail moveTailMoveA(lShort,0);
    
    pimc::moveHead moveHeadMoveB(lShort,1);
    pimc::moveTail moveTailMoveB(lShort,1);
    
    pimc::openMove openA(C, 0, lOpen );
    pimc::closeMove closeA(C, 0, lOpen );
    
    pimc::openMove openB(C, 1, lOpen );
    pimc::closeMove closeB(C, 1, lOpen );

    pimc::createWorm createWormA(C, 0, lShort , 1 );
    pimc::deleteWorm removeWormA(C, 0, lShort , 1);

    pimc::createWorm createWormB(C, 1, lShort , 1 );
    pimc::deleteWorm removeWormB(C, 1, lShort , 1);



    //closeA.setAtomNumberOnClose(N1+N2,{0,1});
    //closeB.setAtomNumberOnClose(N1+N2,{0,1});


    //open.setStartingBead(7);
    //open.setStartingChain(0);

    //close.setStartingBead(7);
    //close.setStartingChain(0);

    //open.setLengthCut(lOpen);
    //close.setLengthCut(lOpen);

    pimc::advanceHead advanceHeadA(lShort,0);
    pimc::recedeHead recedeHeadA(lShort,0);

    pimc::advanceHead advanceHeadB(lShort,1);
    pimc::recedeHead recedeHeadB(lShort,1);
    
    //advanceHead.setFixedLength();
    //recedeHead.setFixedLength();


    pimc::advanceTail advanceTailA(lShort,0);
    pimc::recedeTail recedeTailA(lShort,0);

    pimc::advanceTail advanceTailB(lShort,1);
    pimc::recedeTail recedeTailB(lShort,1);


    pimc::swapMove swapA( lShort , 200 , 0);
    pimc::swapMove swapB( lShort , 200 , 1);

     auto eVEst = std::make_shared<pimc::virialEnergyEstimator>(configurations.nChains() , configurations.nBeads()  );

    auto eEst= std::make_shared< pimc::thermodynamicEnergyEstimator>();


    auto eO= std::make_shared<pimc::scalarObservable>(eEst,std::string("energy") );

    auto eVO= std::make_shared<pimc::scalarObservable>(eVEst,std::string("eV") );


    advanceHeadA.setMaximumParticleNumber( 200 );
    //recedeHeadA.setMinParticleNumber(1);


    advanceHeadB.setMaximumParticleNumber( 200 );
    //recedeHeadB.setMinParticleNumber(1);



    pimc::nConnectedChains nConnectedChains;

    tab.push_back(&levyA,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translateA,0.3,pimc::sector_t::diagonal,"translate");
    tab.push_back(&openA,0.1,pimc::sector_t::diagonal,"open");

    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");

    tab.push_back(&levyA,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translateA,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&closeA,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMoveA,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMoveA,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&advanceHeadA,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHeadA,0.05,pimc::sector_t::offDiagonal,"recedeHead");
    tab.push_back(&swapA,0.1,pimc::sector_t::offDiagonal,"swap");

    tab.push_back(&levyB,0.6,pimc::sector_t::diagonal,"levy");
    tab.push_back(&translateB,0.3,pimc::sector_t::diagonal,"translate");
    tab.push_back(&openB,0.1,pimc::sector_t::diagonal,"open");
    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");


    tab.push_back(&levyB,0.4,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translateB,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&closeB,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMoveB,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMoveB,0.1,pimc::sector_t::offDiagonal,"moveTail");
    tab.push_back(&swapB,0.1,pimc::sector_t::offDiagonal,"swap");
    tab.push_back(&advanceHeadB,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    tab.push_back(&recedeHeadB,0.05,pimc::sector_t::offDiagonal,"recedeHead");


    //configurations.setHeadTail(0,9,-1);


    configurations.fillHeads();

    resetCounters();

    int nTrials=100000;
    int nBlocks=100000;
    
    std::ofstream NOut,l2ShortOut,l2LongOut, ratioOut,particleDistributionOut,wormDistributionOut,lWormOut,nConnectedChainsOut,nRingsOut,magnetizationSquaredOut;


    NOut.open("N.dat");
    l2ShortOut.open("l2Short.dat");
    l2LongOut.open("l2Long.dat");
    particleDistributionOut.open("particleDistribution.dat");
    wormDistributionOut.open("wormDistribution.dat");
    lWormOut.open("lWorm.dat");
    nConnectedChainsOut.open("nConnected.dat");
    nRingsOut.open("nRings.dat");
    ratioOut.open("ratio.dat");

    magnetizationSquaredOut.open("magnetizationSquared.dat");


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

    Real n1Estimator=0;
    Real n2Estimator=0;
    

    Real r=0;

    Real magnetizationSquared=0;

    
    Real nConnectedRingsEstimator=0;
    Real nConnectedChainsEstimator=0;

    for (int i=0;i<nBlocks;i++)
    {

        while (nShort<nTrials and nLong<nTrials and n<nTrials)
        {

            tab.attemptMove(configurations,S,randG);

            if ( configurations.isOpen() )
            {
                r+=1;
            }

            if ( configurations.isOpen(0) )
            {
                
                int l=getWormLength( configurations, 0);

                lWorm+=l;
                l2Short+=accumulateAverageLengthSquare(0,configurations);

                nShort+=1;
                
              /*   if (l < nBeadsWormMax)
                {
                    wormDistribution[l]+=1;
                }

                nConnectedChainsEstimator+=nConnectedChains.count(configurations,0); */

            }
            else if ( not configurations.isOpen() )
            {
                //l2Long+=accumulateAverageLengthSquare( 0,configurations );
                l2Long+=accumulateLengthSquare( configurations , configurations.getGroups()[0].range() , {0,nBeads-1} , geo );

                eVO->accumulate(configurations,S);
                eO->accumulate(configurations,S);

                nLong+=1;
                int currentN1=configurations.nParticles(0);
                int currentN2=configurations.nParticles(1);
                
                n1Estimator+=currentN1;
                n2Estimator+=currentN2;

                magnetizationSquared+=(currentN1 - currentN2) * (currentN1 - currentN2);


                if (currentN1 < nMax )
                {
                    particleDistribution[ currentN1 ] +=1 ;
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
           NOut << i << " " << n1Estimator * 1./nTrials << " " << n2Estimator * 1./nTrials << std::endl << std::flush ;

           magnetizationSquaredOut << i << " " << magnetizationSquared * 1./nTrials  << std::endl << std::flush ;

           for (int iN=0;iN<nMax;iN++)
           {
               particleDistributionOut << i << " " << iN << " " <<  particleDistribution[iN] * 1./nTrials << " " << std::endl << std::flush ;
           }

           nRingsOut << i << " " << nConnectedRingsEstimator/nTrials << std::endl << std::flush ;



           std::fill(particleDistribution.begin(), particleDistribution.end(), 0);

            eO->out(i);
            eO->clear();

            eVO->out(i);
            eVO->clear();

           l2Long=0;
           nLong=0;
           n1Estimator=0;
           n2Estimator=0;
           nConnectedRingsEstimator=0;
           magnetizationSquared=0;


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


#if DIMENSIONS == 3
TEST_F( harmonicTrapTest, caoBernePropagator )
{
    Real C=1e-3;
    int nBeads=10;

    int N=2;
    Real beta=1;
    Real a=0.1;


    std::array<double,DIMENSIONS> lBox = {TRUNCATE_D(1,1,1)};

    SetUp(N,nBeads,beta , lBox );

    //SetUpTwoBodyInteractionHarmonicInTrap();
    //SetUpTwoBodyInteractionGaussian_kernel(10,a);
    SetUpCaoBernePropagator(a);

    //SetUpCaoBernePropagatorTrapped(a);

    //SetUpFreeParticleAction();
    //SetUpTwoBodyInteractionHarmonicInTrap_kernel();

    //SetUpFreeActionWithHardSphereConstraint(a);    
    SetGrandCanonicalEnsamble( 0 );
    SetSeed( time( NULL)  ) ;
    SetRandomMinimumDistance( a, lBox );


    int t0=7;
    int l = int( 0.6* nBeads);
    int lShort=3;
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

// involves levy reconstruction
    pimc::nConnectedChains nConnectedChains;
    tab.push_back(&levy,0.8,pimc::sector_t::diagonal,"levy");
    //tab.push_back(&translate,0.2,pimc::sector_t::diagonal,"translate");
    //tab.push_back(&open,0.1,pimc::sector_t::diagonal,"open");

    //tab.push_back(&createWorm,0.1,pimc::sector_t::diagonal,"createWorm");
    
    tab.push_back(&levy,0.6,pimc::sector_t::offDiagonal,"levy");
    tab.push_back(&translate,0.1,pimc::sector_t::offDiagonal,"translate");
    tab.push_back(&close,0.1,pimc::sector_t::offDiagonal,"close");
    tab.push_back(&moveHeadMove,0.1,pimc::sector_t::offDiagonal,"moveHead");
    tab.push_back(&moveTailMove,0.1,pimc::sector_t::offDiagonal,"moveTail");
    //tab.push_back(&advanceHead,0.05,pimc::sector_t::offDiagonal,"advanceHead");
    //tab.push_back(&recedeHead,0.05,pimc::sector_t::offDiagonal, "recedeHead");
    //tab.push_back(&swap,0.1,pimc::sector_t::offDiagonal,"swap");

    int tTail=4;
    int lWormShort=4;

    //configurations.setHeadTail(0,tTail + lWormShort,tTail -1);


    //configurations.setHeadTail(1,tTail - lCut,-1);
    
    //configurations.join(0,1);
    //configurations.join(1,0);        

    configurations.fillHeads();
    resetCounters();

    int nTrials=1000;
    int nBlocks=100000;

    auto g = std::make_shared<pimc::pairCorrelation>(0,0);
    auto gOb=std::make_shared<pimc::histogramObservable>(g,"pairCorr",1000,0,lBox[0]/2 );
    
    auto eVEst = std::make_shared<pimc::virialEnergyEstimator>(configurations.nChains() , configurations.nBeads()  );

    auto eEst= std::make_shared< pimc::thermodynamicEnergyEstimator>();
    
    auto eO= std::make_shared<pimc::scalarObservable>(eEst,std::string("energy") );

    auto eVO= std::make_shared<pimc::scalarObservable>(eVEst,std::string("eV") );


    int nUnCorrelationSteps = 100;

    configurations.save( "configurationsInitial" , "csv");


    std::ofstream l2Out,dis2Out;
    l2Out.open("l2.dat");
    dis2Out.open("dis2.dat");
    


    Real l2=0;
    Real dis2=0;

    for (int i=0;i<nBlocks;i++)
    {

        int n=0;
        
        while (nClosed<nTrials and n<nTrials)
        {
            for(int n=0;n<nUnCorrelationSteps;n++)
            {
                bool accept = tab.attemptMove(configurations,S,randG);
               

            }
            
            if ( configurations.isOpen() )
            {
               
            }
            else 
            {
                l2+=accumulateAverageLengthSquare(0,configurations);

                dis2+=accumulateAverageParticleDistanceSquared(configurations,{0,M-1},{0,N-1},S.getGeometry());

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


        l2/=nClosed;
        dis2/=nClosed;

        l2Out << i << " " <<  l2 << std::endl;
        dis2Out << i << " " <<  dis2 << std::endl;
    

        nClosed=0;
        l2=0;
        dis2=0;
       }
       

       tab >> std::cout;
       tab.resetCounters();
       configurations.save( "configurations" + std::to_string(i) , "csv");

    }
    
}




#endif



