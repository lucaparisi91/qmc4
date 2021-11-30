
#include "gtest/gtest.h"
#include "geometryPMC.h"
#include "traits.h"
#include "tools.h"
#include <random>
#include <iostream>
#include "../pimc/accelerationStructure.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"
#include "../pimc/hdf5IO.h"
#include "../pimc/toolsPimcTest.h"
#include <filesystem>
#include "../pimc/forces.h"
#include "testConfigurations.h"
#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"

TEST( acceleration, simpleCellNeighbourList )
{  
    int N = 3;
    std::array<Real,DIMENSIONS> lBox = {1,1,1};
    int seed=3456;

    randomGenerator_t randG;

    randG.seed(seed);

    std::uniform_real_distribution<Real> dist(-0.5,0.5);
    std::vector<Real> positions { -0.4, -0.3, 0.25, 0.2, -0.45, 0.49,0.01, 0.32, 0.28};

    std::array<int,DIMENSIONS> nCells {3,3,3};

    acceleration::simpleCellNeighbourList list(nCells,lBox);

    list.add(positions.data(), N ,0, N - 1 );

    list.sanityCheck();

    ASSERT_EQ ( list.cellIndex(0 ) ,15   );
    ASSERT_EQ ( list.cellIndex(1  ) , 18   );
    ASSERT_EQ ( list.cellIndex( 2 ) ,26   );

    ASSERT_EQ(list.nParticles() , 3  );

    for(int iCell=0;iCell<nCells[0]*nCells[1]*nCells[2] ; iCell++ )
    {
        if (iCell==15)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 0 );
        }
        else if (iCell==18)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 1 );
        }
        else if (iCell==26)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 2 );
        }
        else
        {
            ASSERT_EQ(  list.nParticles(iCell) , 0 );
        }       
    }

    positions[2]=0.003;
    positions[2 + N]=-0.1;
    positions[2 + 2* N]=-0.09;

    list.update( positions.data(), N, 2, 2);

    ASSERT_EQ ( list.cellIndex(0 ) ,15   );
    ASSERT_EQ ( list.cellIndex(1  ) , 18   );
    ASSERT_EQ ( list.cellIndex( 2 ) , 13   );

   for(int iCell=0;iCell<nCells[0]*nCells[1]*nCells[2] ; iCell++ )
    {
        if (iCell==15)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 0 );
        }
        else if (iCell==18)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 1 );
        }
        else if (iCell==13)
        {
            ASSERT_EQ(  list.nParticles(iCell) , 1 );
            ASSERT_EQ( list.particle(iCell,0) , 2 );
        }
        else
        {
            ASSERT_EQ( list.nParticles(iCell) , 0 );            
        }   
    }

    int nParticlesNeighBours=0;
    int iCell=13;


    for(int i=0;i < list.nCellNeighbours( iCell ) ;i++)
    {
        auto iCellNext=list.neighbourCell(iCell,i);

        nParticlesNeighBours+=list.nParticles(iCellNext);
    }

    ASSERT_EQ( nParticlesNeighBours, 3   );

    iCell=18;
    nParticlesNeighBours=0;   

    ASSERT_EQ( list.nCellNeighbours(iCell) , 27  );


    for(int i=0;i < list.nCellNeighbours( iCell ) ;i++)
    {
        auto iCellNext=list.neighbourCell(iCell,i);
        nParticlesNeighBours+=list.nParticles(iCellNext);
    }
    ASSERT_EQ( nParticlesNeighBours, 3   );


}


TEST_F( configurationsTest, acceleration )
{  
    int N=100;
    int M=100;
    std::array< Real, DIMENSIONS> lBox {1,1,1};
    std::array< int, DIMENSIONS> nCells {10,10,10};

    SetUp( {N},M,1.0);

    SetRandom(lBox);

    pimc::cellNeighbourList acc;
    acc.setGeometry(geo);
    acc.setNBeads(M);
    acc.setNCells(nCells);

    acc.build();
    
    acc.add(configurations.dataTensor(),{0,M}, {0,N-1} );


}

