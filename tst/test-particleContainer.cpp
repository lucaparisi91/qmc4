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

#include "testConfigurations.h"


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
        cell.push( i)  ;
    }

    ASSERT_EQ(cell.size() , N );

    const auto & indices=cell.indices();

    for(int i=0;i<N;i++)
    {
        for(int d=0;d<getDimensions();d++)
        {
            //ASSERT_NEAR(cellPositions(i,d) , positions(i,d) , TOL);
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
            //ASSERT_NEAR(cellPositions(i,d) , positions(indices[i],d) , TOL);
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

auto countDistances( const Eigen::Tensor<Real,3> & data,Real cutOff , int t , int iParticle , const pimc::range_t & particleRange , const pimc::geometry_t & geo)
{
    Real cutOffSquared=cutOff*cutOff;

    size_t n=0;
    Real diff=0;
    for(int i=particleRange[0];i<=particleRange[1];i++)
    {
        if (i==iParticle)
        {
            continue;
        }
        Real r2=0;
        Real r2Next=0;

        for(int d=0;d<getDimensions();d++)
        {
            Real diff=geo.difference(data(i,d,t) - data(iParticle,d,t),d);
            Real diffNext=geo.difference(data(i,d,t+1) - data(iParticle,d,t+1),d);

            r2+=diff*diff;
            r2Next+=diffNext*diffNext;
        }

        if (not ( (  r2>=cutOffSquared  ) or   (r2Next>=cutOffSquared  ) ) )
        {
            n+=1;
            //std::cout << t << " " << i << " " << iParticle << std::endl;
        }
       
    }

    return n;

}




TEST( particleDistances , basic)
{
    int N=500;
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
        //std::cout << t << " ";
        for(int i=0;i<nFilteredParticles[t];i++)
        {
            //std::cout << filteredParticles(i,t) << " ";
        }
        //std::cout << std::endl;
    }



    pimc::twoBodyDistancesBuffer distanceBuffer(N,nBeads);
    pimc::geometry_t geo(lBox[0],lBox[1],lBox[2]);
    

    distanceBuffer.buildDistanceList( positions,griddedParticles,geo,{0,nBeads-1},0);

    //std::cout  << "...distancesBuffer.." << std::endl;

    

    for (int t=0;t<nBeads;t++)
    {
        auto n1=distanceBuffer.getNDistances()[t];
        auto n2 = countDistances( positions , griddedParticles.cellLength(0),t,iParticle,{0,N-1},geo  );

        ASSERT_EQ(n1,n2);


    }




     /*
    TODO


    Some check on filtered distances

    */




}

#include "../pimc/pairProductKernel.h"
#include "../pimc/actionTwoBodyConstructor.h"

class caoBerneGridded : public ::testing::Test
{
protected:

    void SetUp()
    {
        nCells={3,3,3};
        lBox= { 1,1,1};
        beta=1;
        nBeads=10;
        N=10;
        scatteringLength=0.01;
        _doCutOff=false;
        nBurnsIn=0;


    };

    void SetBox( std::array<Real,getDimensions()> l)
    {
        lBox=l;
    };

    void SetBeta(Real beta_)
    {
        beta=beta_;
    }


    void SetNCells( std::array<size_t,getDimensions()> nCells_)
    {
        nCells=nCells_;
    };

    void SetScatteringLength(Real a)
    {
        scatteringLength=a;

    }

    void SetN( int N2)
    {
        N=N2;
    };

    void SetNBeads( int M)
    {
        nBeads=M;
    }


    void SetNBurnsIn( int M)
    {
        nBurnsIn=M;
    }


    void run()
    {

        Real timeStep=beta/nBeads;


        randomGenerator_t randG( time(NULL));



        pimc::pimcConfigurations configurations(nBeads,getDimensions(), { {0,N-1,N+1,1.0} } );

        pimc::geometryPBC_PIMC geo({lBox[0],lBox[1],lBox[2] });





        int seed=567;

        configurations.setRandom( lBox , randG );

        configurations.fillHeads();

        auto sTwoBodyCreator = std::make_shared<pimc::actionTwoBodyConstructor >() ;
        sTwoBodyCreator->setTimeStep(timeStep);
        sTwoBodyCreator->setGeometry(geo);
        sTwoBodyCreator->setNBeads(nBeads);

        sTwoBodyCreator->setNMaxParticles(N);

        sTwoBodyCreator->registerGreenFunction<pimc::testingCaoBernePropagator>("caoBerne");



        auto j = R"(
        {
                "kind": "twoBody",
                "groupA": 0,
                "groupB": 0,
                "greenFunction": {
                    "a": 0.01,
                    "kind": "caoBerne"
                }
            }

                )"_json;

        j["greenFunction"]["a"]=scatteringLength;



        Real cutOff=lBox[0]*2;
        if (_doCutOff)
            {
                cutOff=lBox[0]*1./nCells[0];
            }
    

       
        j["greenFunction"]["cutOff"]=cutOff;


        auto sG = sTwoBodyCreator->create( j);

        while (not sG->checkConstraints(configurations) )
        {
            configurations.setRandom( lBox , randG );
        }

        
        auto sT= std::make_shared<pimc::kineticAction>(timeStep, configurations.nChains() , nBeads  , geo);

        auto S = std::make_shared<pimc::firstOrderAction>(sT, sG);

        int l= int(nBeads*0.6);
        pimc::levyMove levyMove(l,0);

        
         for(int n=0;n<nBurnsIn;n++)
        {
            levyMove.attemptMove(configurations,*S,randG);
        }




        pimc::range_t timeRange {0,nBeads-1};
        int iParticle=0;



        auto G = std::make_shared<pimc::caoBernePropagator>(timeStep,scatteringLength);

        auto meshKernel=std::make_shared<pimc::pairProductMeshKernel<pimc::caoBernePropagator> >(G,N,nBeads);



        meshKernel->setGeometry(geo);
        Real a = j["greenFunction"]["a"].get<Real>();

        p=0;
        p2=0;
        int nSteps=1000;

        for(int n=0;n<nSteps;n++)
        {
            for(int nn=0;nn<10;nn++)
            {
                levyMove.attemptMove(configurations,*S,randG);
            }
            
            pimc::linkedCellParticles linkedCellList( nCells, lBox);    

           
 
            linkedCellList.setCapacity(N,nBeads);
            linkedCellList.add( configurations.dataTensor(), {0,nBeads-1} , {0,N-1} );

            p += sG->evaluate(configurations, timeRange, iParticle);
        
            p2+=meshKernel->evaluate( configurations.dataTensor() , linkedCellList,{0,nBeads} , iParticle );

        }
        










         



    }



    auto getNonGriddedResult()  {return p;}
    auto getGriddedResult() {return p2;}

    void SetDoCutOff( bool doCutOff)  {_doCutOff=doCutOff;} 
    
    private:

    std::array<size_t,getDimensions()> nCells;
    std::array<Real,getDimensions()> lBox;
    Real beta;
    int nBeads;
    int N;
    Real scatteringLength;

    Real p;
    Real p2;

    bool _doCutOff;
    int nBurnsIn;

};


TEST_F( caoBerneGridded , nCells3)
{

    Real na3=1e-6;
    std::vector<int> Ns { 100 , 200 , 300 , 500 };
    
    Real T=0.1;
    Real zeta=2.612;

    SetBeta(1);


    for (auto N : Ns )
    {
        Real V = N/zeta * std::pow(T,3/2.);

        Real l=std::pow(V,1/3.);
        SetScatteringLength( std::pow(na3/(N/V)  ,1./3) );

        SetBox( {l,l,l});
        SetN(N);

        run();

        auto p = getNonGriddedResult();
        auto p2=getGriddedResult();

        ASSERT_EQ(p,p2);

    }

}
TEST_F( caoBerneGridded , cutOff)
{
    Real na3=1e-6;
    std::vector<int> Ns { 100  };
    std::vector<int> nCells { 3 , 4 , 5 , 6  };

    Real T=0.1;
    Real zeta=2.612;

    Real beta=1./(2*M_PI);

    int nBeads=100;


    SetNBeads(nBeads);
    SetBeta(beta);

    SetDoCutOff(true);




    for (auto N : Ns )
    {
        Real V = N/(zeta * std::pow(T,3/2.) );

        Real l=std::pow(V,1/3.);
        
        Real a=std::pow(na3/(N/V),1./3 )  ;
        SetScatteringLength( a  ) ;


        for( auto nCell : nCells)
        {

        
            
            SetBox( {l,l,l});
            SetN(N);
            SetNCells( {nCell,nCell,nCell});

            run();

            auto p = getNonGriddedResult();
            auto p2=getGriddedResult();


            ASSERT_NEAR(p,p2, 1e-9);



        }

    }

}

TEST_F( caoBerneGridded , accuracy)
{
    Real na3=1e-6;
    std::vector<int> Ns { 100  };
    std::vector<int> nCells { 3 , 4 , 5 , 6 , 8 , 10 };


    Real T=0.1;
    Real zeta=2.612;

    Real beta=1./(2*M_PI);

    int nBeads=100;

    SetNBeads(nBeads);
    SetBeta(beta);

    SetDoCutOff(false);
    
    
    SetNBurnsIn(1000);

    for (auto N : Ns )
    {
        Real V = N/(zeta * std::pow(T,3/2.) );

        Real l=std::pow(V,1/3.);
        
        Real a=std::pow(na3/(N/V),1./3 )  ;
        SetScatteringLength( a  ) ;

       
        std::cout << "eps : " <<  ( sqrt(beta/nBeads) + a)/l << std::endl;
        std::cout << "l : " <<  l << std::endl;
        std::cout << "tau : " <<  beta/nBeads << std::endl;
        std::cout << "a : " <<  a << std::endl;


        for( auto nCell : nCells)
        {
        
            
            SetBox( {l,l,l});
            SetN(N);
            SetNCells( {nCell,nCell,nCell});
            
            Real p=0;
            Real p2=0;

           
            run();
            p= getNonGriddedResult();
            p2=getGriddedResult();

            


            std::cout << nCell << " " <<  p << " "<< p2 << std::endl;


        }

    }


}




