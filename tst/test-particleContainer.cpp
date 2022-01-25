#include "integrationTests.h"
#include "../pimc/actionTwoBody.h"
#include "../pimc/potentialKernel.h"
#include "../pimc/pairProductKernel.h"
#include "../pimc/propagators.h"
#include "testConfigurations.h"
#include "../pimc/pimcObservables.h"
#include <filesystem>
#include "../pimc/particleContainer/particleContainer.h"
#include "../pimc/particleContainer/particleDistancesList.h"
#include "../pimc/pairProductMeshKernel.h"



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


TEST( tools , listIntersection)
{
    std::vector<int> L1 { 1 , 3 , 6 ,8 };
    std::vector<int> L2 { 1, 2 , 6  , 7, 8 , 9  };

    int n = pimc::intersect( L1.begin(),L1.end(),L2.begin(),L2.end() );

    ASSERT_EQ( n , 3);

    ASSERT_EQ(L1[0],1);
    ASSERT_EQ(L1[1],6);
    ASSERT_EQ(L1[2],8);

}


TEST( particleDistances , basic)
{
    int N=10;
    int nBeads=10;


    pimc::cell cell{};
    cell.setCapacity( N );

    


    std::array<Real,getDimensions() > lBox { 1 , 1 , 1};

    Eigen::Tensor<Real , 3 > positions( N, getDimensions() , nBeads +1 );

    std::default_random_engine randG;
    std::uniform_real_distribution<double> distribution(-0.5,0.5);

    for( int t=0;t<nBeads;t++)
        for(int i=0;i<N;i++)
        {
            for(int d=0;d<getDimensions(); d++)
            {
                positions(i,d,t)=distribution(randG)*lBox[d] ;
            }
        }

    for(int i=0;i<N;i++)
        {
            for(int d=0;d<getDimensions(); d++)
            {
                positions(i,d,nBeads)=positions(i,d,0) ;
            }
        }

    std::array<size_t,getDimensions()> nCells { 4,4,4};

    pimc::linkedCellParticles griddedParticles(nCells,lBox);

    griddedParticles.setCapacity(N,nBeads);
    
    griddedParticles.add( positions, {0, nBeads -1} , {0,N-1} );

    pimc::twoBodyPairsList pairList(N,nBeads);
    pairList.buildParticleList(griddedParticles, 0, {0,nBeads-1});



    const auto & filteredParticles = pairList.getPairList();
    const auto & nFilteredParticles = pairList.sizes();
    /*
    TODO

    Some check on particles filtered to implement
    */

    int iParticle=0;

    for(int t=0;t<nBeads;t++)
    {

        for(int i=0;i<nFilteredParticles[t];i++)
        {
            std::cout << filteredParticles(i,t) << " ";
        }
        std::cout << std::endl;
    }



    pimc::twoBodyDistancesBuffer distanceBuffer(N,nBeads);
    pimc::geometry_t geo(lBox[0],lBox[1],lBox[2]);
    

    distanceBuffer.buildDistanceList( positions,griddedParticles,geo,{0,nBeads-1},0);

     /*
    TODO


    Some check on filtered distances

    */

   
}

