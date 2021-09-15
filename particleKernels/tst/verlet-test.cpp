#include "gtest/gtest.h"
#include "cellNeighbourList.h"
#include <random>
#include <cmath>
#include "timers.h"
#include "testUtils.h"


using namespace particleKernels;


TEST(metaprogramming,explore)
{
    
};




Real differencePBC(Real t, Real lBox, Real lBoxInverse ) {return ( t - std::round(t*lBoxInverse )*lBox);}


Real difference(Real t, Real lBox, Real lBoxInverse ) {return  t ;}



template<int D,class V_t>
double twoBodyPotential(double * particles,int iStart, int iEnd, int N,const V_t & V,double cutOff,double * lBox)
{
    double sumCheck=0;
    for (int i=iStart;i<=iEnd;i++)
    {
        for (int j=0;j<i;j++)
        {
            double dis=0;

            for (int d=0;d<D;d++)
            {
                double tmp=differencePBC((particles[i+N*d] - particles[j+N*d]),lBox[d],1./lBox[d]);

                dis+=tmp*tmp;
            }          

            dis=std::sqrt(dis);

            if (dis <= cutOff)
            {
                sumCheck+=V(dis);
            }

        }

    }

    return sumCheck;
}










TEST(cellNneighbourList,configure)
{   
    const int N=100000;

    constexpr int D = 3;

    int iStart=0;
    int iEnd=N-1;

    const int nTrials = 10000;


    Real density = 1e-3;

    Real l = std::pow((N/density) , 1./3 );

    Real cutOff=1;
    int nCellMax=100;




    ASSERT_GT(l,cutOff);


    std::vector<Real> lBox{l,l,l};
    std::vector<int> nCells{1,1,1};

    

    for (int d=0;d<D;d++)
    {
        nCells[d]=std::min(  (int)( lBox[d]/(cutOff))  , nCellMax);
        ASSERT_GT(nCells[d] , 0 );

    }

    

    // generates particle configurations

    std::allocator<Real> alloc;

    double * particles = alloc.allocate(N*D);

    std::default_random_engine generator;
    std::uniform_real_distribution<double> unifRandom(0.0,1.0);


   generateRandomParticlesInBox(particles,iStart,iEnd,N ,D,  generator, lBox.data() );


    gaussianPotential V(0.1);




    ADD_TIMER("listConstructor");
    ADD_TIMER("listTwoBodyPotential");
    ADD_TIMER("twoBodyPotential");
    ADD_TIMER("listInitialization");

    START_TIMER("listInitialization");
    simpleCellNeighbourList accStructure(nCells,lBox);
    STOP_TIMER("listInitialization");


    accStructure.reportCellNeighbourTopology() ;

    accStructure.sanityCheck();

    accStructure.list(particles,N,iStart,iEnd);

    double sumTot=0;
    double sumTotCheck=0;


    for (int i=0;i<nTrials;i++)
    {

        generateRandomParticlesInBox(particles,iStart,iEnd,N ,D,  generator, lBox.data() );

#ifndef NDEBUG
        START_TIMER("twoBodyPotential");
        double sumCheck=twoBodyPotential<3,decltype(V)>(particles,iStart,iEnd,N,V,cutOff,lBox.data() );
        STOP_TIMER("twoBodyPotential");
#endif



        START_TIMER("listConstructor");
        accStructure.updateList(particles,N,iStart,iEnd);
        STOP_TIMER("listConstructor");
        #ifndef NDEBUG
        accStructure.reportParticleDistribution(particles,0,N-1,N) ;
        #endif


        START_TIMER("listTwoBodyPotential");
        Real sum=accStructure.twoBodyPotential(particles,iStart,iEnd,N,V,cutOff);
        STOP_TIMER("listTwoBodyPotential");
        
        sumTot+=sum;
#ifndef NDEBUG
        sumTotCheck+=sumCheck;
#endif
        

    }

    std::cout << "SumAcc: " << sumTot << std::endl;
#ifndef NDEBUG
    std::cout << "SumCheck: " << sumTotCheck << std::endl;

    ASSERT_NEAR(sumTot,sumTotCheck,1e-7);
#endif


    

    
    std::cout << "nCells : "<< nCells[0] << std::endl;


    
    std::cout << timers::getInstance().report() << std::endl;


#ifndef NDEBUG


    std::cout << "SpeedUp : " <<
    (timers::getInstance()["twoBodyPotential"].timeElapsed()) /
    ( 
        timers::getInstance()["listConstructor"].timeElapsed() + timers::getInstance()["listTwoBodyPotential"].timeElapsed() 
    ) << std::endl;
#endif

    alloc.deallocate(particles,N*D);
}






