#include "integrationTests.h"
#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"
#include "testConfigurations.h"
#include "../pimc/pimcObservables.h"
#include <filesystem>
#include "../pimc/particleContainer/particleContainer.h"

namespace fs = std::filesystem;


TEST( cell , basic)
{
    int N=10;
    pimc::cell cell{};
    cell.setCapacity( N );
    
    std::array<Real,getDimensions() > lBox { 1 , 1 , 1};

    Eigen::Tensor<Real , 2 > positions( N, getDimensions() );

    std::default_random_engine randG;
    std::uniform_real_distribution<double> distribution(-0.5,0.5);

    for(int i=0;i<N;i++)
    {
        for(int d=0;d<N;d++)
        {
            positions(i,d)=distribution(randG)*lBox[d] ;
        }
    }
    

    cell.setCapacity(N);

    ASSERT_EQ(cell.capacity() , N);

    ASSERT_EQ(cell.size() , 0 );

    for(int i=0;i<N;i++)
    {
        cell.push( i,positions(i,0) ,positions(i,1),positions(i,2)  );
    }

    ASSERT_EQ(cell.size() , N );

    const auto & cellPositions=cell.positions();
    const auto & indices=cell.indices();

    for(int i=0;i<N;i++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            ASSERT_NEAR(cellPositions(i,d) , positions(i,d) , TOL);
            ASSERT_EQ(indices[i],i);
        }
    }

    int i=1;
    int j=2;

    cell.swap(i,j);
    ASSERT_EQ( indices[j],i);
    ASSERT_EQ( indices[i],j);


    

    for(int i=0;i<N;i++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            ASSERT_NEAR(cellPositions(i,d) , positions(indices[i],d) , TOL);
        }
    }

    cell.pop();

    ASSERT_EQ(cell.size(), N-1);
    ASSERT_EQ(cell.capacity(), N);


    ASSERT_EQ(cell.neighbours().size() , 26);
    ASSERT_EQ(cell.displacements().dimensions()[0] , 26);
    ASSERT_EQ(cell.displacements().dimensions()[1] , getDimensions());

}

TEST( simpleCellNeighbourList , basic)
{
    int N=10;
    pimc::cell cell{};
    cell.setCapacity( N );
    
    std::array<Real,getDimensions() > lBox { 1 , 1 , 1};

    Eigen::Tensor<Real , 2 > positions( N, getDimensions() );

    std::default_random_engine randG;
    std::uniform_real_distribution<double> distribution(-0.5,0.5);

    for(int i=0;i<N;i++)
    {
        for(int d=0;d<getDimensions(); d++)
        {
            positions(i,d)=distribution(randG)*lBox[d] ;
        }

    }

    std::array<size_t,getDimensions()> nCells { 10,10,10};

    pimc::simpleCellNeighbourList linkedCells( nCells, lBox);

    ASSERT_EQ(linkedCells.size(),0);


    linkedCells.setCapacity(N);

    std::array<Real,3> testPos {0.105, 0.25,0.48};
    linkedCells.add( 0, testPos[0] , testPos[1], testPos[2] );

    ASSERT_EQ(linkedCells.size() , 1 );

    ASSERT_EQ(linkedCells.getCells().size(),nCells[0]*nCells[1]*nCells[2]  );

    std::array<int,3> testCellIndex {6,7,9};


    const auto & updatedCell=linkedCells.getCell( testCellIndex[0] , testCellIndex[1],testCellIndex[2] );

    ASSERT_EQ( updatedCell.size(), 1 );
    ASSERT_EQ(updatedCell.getCellIndex()[0],testCellIndex[0]);
    ASSERT_EQ(updatedCell.getCellIndex()[1],testCellIndex[1]);
    ASSERT_EQ(updatedCell.getCellIndex()[2],testCellIndex[2]);

    ASSERT_EQ(updatedCell.indices()[0],0);



    const auto & positionsCell = updatedCell.positions();

    ASSERT_EQ( positionsCell(0,0), testPos[0] );
    ASSERT_EQ( positionsCell(0,1), testPos[1] );
    ASSERT_EQ( positionsCell(0,2), testPos[2] );



    linkedCells.remove(0);
    ASSERT_EQ(updatedCell.size(),0);

    const auto & neighbours=updatedCell.neighbours();

    for (int i=0;i<neighbours.size();i++)
    {

        const auto & neighbour = neighbours[i];

        ASSERT_FALSE(  neighbour == nullptr);

        for ( int d=0;d<getDimensions();d++)
        {
            
            //std::cout << neighbour->getCellIndex()[d] << " ";

            Real dis=neighbour->getCellIndex()[d] - testCellIndex[d];
            Real delta = nCells[d]/2.0;

            Real displacement=0;

            if (dis>= delta )
            {
                dis-=nCells[d];
                displacement=lBox[d];
            }

            if(dis<= -delta )
            {
                dis+=nCells[d];
                displacement=-lBox[d];
            }

            ASSERT_LE(  std::abs( dis ) , 1   );
            ASSERT_NEAR( updatedCell.displacements()(i,d) , displacement, TOL    );

        }
        std::cout<<std::endl;
    }

}

