#include "testConfigurations.h"


void configurationsTest::SetRandom(const std::array<Real,DIMENSIONS> & lBox )
    {
        std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

        auto & data=configurations.dataTensor();
        for (int t=0;t<data.dimensions()[2];t++)
            for (int i=0;i<data.dimensions()[0];i++)
                for  (int d=0;d<DIMENSIONS;d++)
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

#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"
#include "../pimc/actionTwoBody.h"

void configurationsTest::SetUpTwoBodyInteractionGaussian_pairProductKernel( Real V0 , Real alpha )
    {
        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);
    
        auto V = pimc::makeIsotropicPotentialFunctor(
         [V0,alpha](Real r) {return V0*exp(-alpha*r*r) ;} ,
         [V0,alpha](Real r) {return -V0*2*alpha*r*exp(-alpha*r*r)  ;} );

        using propagator_t=pimc::primitivePropagator<decltype(V)>;


        auto G= std::make_shared< propagator_t >(timeStep,V);
        auto pairKernel = std::make_shared< pimc::pairProductKernel<propagator_t> >( G );
        pairKernel->setTimeStep(timeStep);
        pairKernel->setGeometry(geo);
        auto SPP = std::make_shared<pimc::actionTwoBody>();
        SPP->setTimeStep(timeStep);
        SPP->setGeometry(geo);
        SPP->setKernel(pairKernel);
        SPP->setSets({0,0});


       S= pimc::firstOrderAction(sT,  SPP);

    }







void configurationsTest::SetUpTwoBodyInteractionHarmonicInTrap()
{

    SetUpTwoBodyInteractionHarmonicInTrap( { {0,0} });

}


 void configurationsTest::SetUpTwoBodyInteractionHarmonicInTrap( const std::vector<std::pair<int,int> > & sets)
    {
        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);

        int order=2;

        auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );
        
        std::vector<std::shared_ptr<pimc::action> > Vs;

        for ( int set=0;set<configurations.getGroups().size() ; set++)
        {
            auto  sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,set,order);
            Vs.push_back(sOneBody);

        }
        

        for ( auto [ setA, setB ] : sets )
        {
            auto sV2B=std::make_shared<pimc::potentialActionTwoBody<decltype(V)>  >(timeStep,configurations.nChains(),configurations.nBeads(),V,geo,setA,setB);

            Vs.push_back(sV2B);

        }

        std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

        S= pimc::firstOrderAction(sT,  sV);

    }

    void configurationsTest::SetUpNonInteractingHarmonicAction(int order)
    {
         std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , M  , geo);



     auto V = pimc::makeIsotropicPotentialFunctor(
         [](Real r) {return 0.5*(r*r) ;} ,
         [](Real r) {return r  ;} );

    
    std::vector<std::shared_ptr<pimc::action> > Vs;

    for( int set=0;set<configurations.getGroups().size() ;set++  )
    {
        std::shared_ptr<pimc::action> sOneBody=std::make_shared<pimc::potentialActionOneBody<decltype(V)> >(timeStep,V ,geo,set,order);
        Vs.push_back( sOneBody);
    }


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
                    Real diffd= data( i ,d, t+1) - data(i,d,t) ;
                    l2+=diffd*diffd;
                }
            }
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


Real accumulateAverageParticleDistanceSquared( const  pimc::configurations_t & configurations , const std::array<int,2> & timeRange , const std::array<int,2> & particleRange , const pimc::action::geometry_t & geo)
{
    const auto & data = configurations.dataTensor();
    Real sum=0;

    for(int t=timeRange[0];t<=timeRange[1];t++){ 
        for(int i=particleRange[0];i<=particleRange[1];i++) { 
            for(int j=particleRange[0];j<i;j++) {
                for(int d=0;d<DIMENSIONS;d++)
                {
                    Real diff = geo.difference( data(i,d,t) - data(j,d,t) , d );
                    sum+=diff*diff;
                }
            }
        }
    }

    auto M = timeRange[1] - timeRange[0] + 1;
    auto N= particleRange[1] - particleRange[0] + 1;

    return sum/(N*M);
}







std::array<Real,getDimensions()> meanBeadFixedLengths(int iChainBegin , int iChainEnd,  int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
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

std::array<Real,getDimensions()> varianceBeadFixedLengths(int iChain , int t0, int t1, int i, Real timeStep, const pimc::configurations_t & configurations)
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



