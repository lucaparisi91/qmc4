#include "gtest/gtest.h"
#include <random>
#include <cmath>
#include "lattice.h"
#include <cmath>
#include <random>
#include "particleData.h"
#include "twoBodyDistances.hpp"


#include "utils.h"








using namespace particleKernels;

TEST(linkedCell,lattice)
{

    Real cut_off=1;
    constexpr int D = 3;
    size_t N = 10000;
    Real density = 1e-6;

    using testParticleData_t = particleData3D;

    

    Real l=std::pow(N/density,1./3);


    std::array<Real, D> lBox={l,l,l};
    std::array<size_t,D> nBoxes;
    std::array<Real,D> leftEdge={-lBox[0]/2.,-lBox[1]/2.,-lBox[2]/2.};
    std::array<Real,D> rightEdge={lBox[0]/2.,lBox[1]/2.,lBox[2]/2.};

    for (int d=0;d<D;d++)
    {
        nBoxes[d]=std::min(std::pow(N,1./3),std::floor(lBox[d]/(cut_off) ) );
    }



    auto lattptr= std::shared_ptr<lattice>( new  lattice(nBoxes, leftEdge,rightEdge));


    //ASSERT_EQ(lattptr->lowIndex(0),1);
    //ASSERT_EQ(lattptr->highIndex(0),10);

    lattptr->checkNeighbourIndexing();


    double * positions = (Real *)aligned_alloc( 128, N*D*sizeof(Real)  );
    double * positionsPBC = (Real *)aligned_alloc( 128, N*D*sizeof(Real)  );


    std::default_random_engine generator(12);


    utils::initRandom(positions,3,0,N-1,N,generator,-5,5);


    utils::restrictToBox(positionsPBC,positions,0,N-1,N,leftEdge,lBox);


    utils::gaussianInteractionWithCutOff gauss(2,1);
    utils::tetaInteraction teta(1,1);


    std::cout << "Direct two body differences" << std::endl;
    Real vCheck = utils::twoBodyDistancesIsotropicReduction<utils::gaussianInteractionWithCutOff,3>(positions,gauss,N,lBox);


    std::cout << "Sum: " << vCheck << std::endl;

    std::cout << "Constructing data structure" << std::endl;
    
    testParticleData_t particleAcc(lattptr);

    std::cout << "Listing particles" << std::endl;
    particleAcc.list(positionsPBC,0,N-1,N);

    particleAcc.checkPBCGhostCells();

    std::cout << "Two body distances accellerated" << std::endl;


    auto v = twoBodyDistancesIsotropicReduction(particleAcc,gauss,positionsPBC,0,N-1,N);

    ASSERT_NEAR(v,vCheck,1e-6);

    std::uniform_int_distribution<int> sampleParticle(0,N-1);

    for (int n=0;n<10;n++)
    {
        int iA=sampleParticle(generator);
        int iB=sampleParticle(generator);

        //int iStart=std::min(iA,iB);
        //int iEnd=std::max(iA,iB);

        int iStart=0;
        int iEnd=N-1;

        std::cout << "Initialize" << std::endl;
        utils::initRandom(positions,3,iStart,iEnd,N,generator,-5,5);

         std::cout << "Restrict to box" << std::endl;
        utils::restrictToBox(positionsPBC,positions,iStart,iEnd,N,leftEdge,lBox);

        std::cout << "Direct two body differences" << std::endl;
        Real vCheckUpdate = utils::twoBodyDistancesIsotropicReduction<utils::gaussianInteractionWithCutOff,3>(positions,gauss,N,lBox);

        std::cout << "updating list" << std::endl;
        particleAcc.updateList(positionsPBC,iStart,iEnd,N);

        std::cout << "Two body differences acc." << std::endl;
        Real vUpdate=twoBodyDistancesIsotropicReduction(particleAcc,gauss,positionsPBC,iStart,iEnd,N);

        ASSERT_NEAR(vUpdate,vCheckUpdate,1e-6);

    }

    free(positions);
    free(positionsPBC);

};
