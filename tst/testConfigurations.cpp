#include "testConfigurations.h"

void configurationsTest::SetRandom(const std::array<Real,DIMENSIONS> & lBox )
    {
        std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

        auto & data=configurations.dataTensor();
        for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<getDimensions();d++)
            {
                data(i,d,t)=(uniformDistribution(randG)-0.5 )*lBox[d];
            }
        configurations.fillHeads();
    
    }

 void configurationsTest::SetUp( int N_ , int M_ , Real Beta_,std::array<Real,getDimensions()> lBox) {
        N=N_;
        M=M_;
        Beta=Beta_;

        int seed= 356;
        int buffer=2000;

        int nChains=N + buffer;

        timeStep=Beta/M;
        pimc::particleGroup groupA{ 0 , N-1, nChains -1 , 1.0};


        configurations=pimc::pimcConfigurations(M , getDimensions() , {groupA});

        geo = pimc::geometryPBC_PIMC(lBox[0],lBox[1],lBox[2]);

        randG=randomGenerator_t(seed);

        SetRandom();

    }



void configurationsTest::SetUp( std::vector<int> Ns , int M_ , Real Beta_,std::array<Real,getDimensions()> lBox) {
        N=std::accumulate(Ns.begin(),Ns.end(),0);
        M=M_;
        Beta=Beta_;

        int seed= 356;
        int buffer=1000;

        int nChains=N + buffer*Ns.size();

        timeStep=Beta/M;
        std::vector<pimc::particleGroup> groups;

        int i=0;
        for (int iSet=0;iSet<Ns.size();iSet++)
        {
            pimc::particleGroup group{ i , i + Ns[iSet]-1, i + Ns[iSet] -1 + buffer , 1.0};
            groups.push_back(group);
            i+= Ns[iSet] + buffer;
        }
        

        configurations=pimc::pimcConfigurations(M , getDimensions() , groups);

        geo = pimc::geometryPBC_PIMC(lBox[0],lBox[1],lBox[2]);

        randG=randomGenerator_t(seed);

        SetRandom();

    }



void configurationsTest::SetUpFreeParticleAction()
    {
         std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

    std::shared_ptr<pimc::action> sV= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);


    
    S= pimc::firstOrderAction(sT,  sV);

    }


void configurationsTest::SetUpTwoBodyInteractionHarmonic()
    {
        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

        

        auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );

        auto sV=std::make_shared<pimc::potentialActionTwoBody<decltype(V)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V,geo,0,0);

       S= pimc::firstOrderAction(sT,  sV);

       
    }


void configurationsTest::SetUpTwoBodyInteractionGaussian( Real V0 , Real alpha )
    {
        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);


        auto V = pimc::makeIsotropicPotentialFunctor(
         [V0,alpha](Real r) {return V0*exp(-alpha*r*r) ;} ,
         [V0,alpha](Real r) {return -V0*2*alpha*r*exp(-alpha*r*r)  ;} );

        auto sV=std::make_shared<pimc::potentialActionTwoBody<decltype(V)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V,geo,0,0);

       S= pimc::firstOrderAction(sT,  sV);

    }

 void configurationsTest::SetUpTwoBodyInteractionHarmonicInTrap()
    {
        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

        int order=2;

        auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );

        auto sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V,geo,0,0);

        auto  sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,order);

        std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody,sV2B};

        std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

        S= pimc::firstOrderAction(sT,  sV);

    }

    void configurationsTest::SetUpNonInteractingHarmonicAction(int order)
    {
         std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);




     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );



    std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,order);
    
    
    std::vector<std::shared_ptr<pimc::action> > Vs = {sOneBody};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    S= pimc::firstOrderAction(sT,  sV);

    }



Real accumulateLengthSquare(const pimc::configurations_t & configurations, const std::array<int,2> & particleRange, std::array<int,2> timeRange,const geometryPBC & geo)
    {
        const auto & data = configurations.dataTensor();
        Real l2=0;

        for(int i=particleRange[0]; i<=particleRange[1];i++)
        {
            const auto & chain = configurations.getChain(i);

            for (int t=std::max(chain.tail + 1, timeRange[0]);t<=std::min(chain.head , timeRange[1] ) ;t++)
            {
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diffd=geo.difference( data( i ,d, t+1) - data(i,d,t) ,d );
                    l2+=diffd*diffd;
                }
            }
        }

        return l2;
    }

