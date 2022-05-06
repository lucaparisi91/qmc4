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
#include "../pimc/actionTwoBodyMesh.h"

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
            Real diffNext= -(  - diff - ( data(i,d,t+1) - data(i,d,t) ) +  data(iParticle,d,t+1) - data(iParticle,d,t)); 

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


TEST( gridded_evaluation, caoBerne )
{

}



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
        _doCutOff=true;
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
    auto getNBeads() const {return nBeads; }





    void SetNCells( std::array<size_t,getDimensions()> nCells_)
    {
        nCells=nCells_;
    };
    auto getNCells( )
    {
        return nCells;
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

   
    

    auto getTimeStep()
    {
        return beta/nBeads;
    }

    auto createGeometry()
    {
        return pimc::geometryPBC_PIMC({lBox[0],lBox[1],lBox[2] });

    }


    

    auto createLinkedCellParticles()
    {
        auto linkedCellList=std::make_shared<pimc::linkedCellParticles> ( nCells, lBox);    
        linkedCellList->setCapacity(N,nBeads);
        return linkedCellList;
    }


    auto createConfigurations( randomGenerator_t & randG)
    {
        pimc::pimcConfigurations confs(nBeads,getDimensions(),{{0,N-1,N+1,1.0} });

        confs.setEnsamble(pimc::ensamble_t::grandCanonical);

        generateRandomMinimumDistance(  confs, scatteringLength,randG, createGeometry());


        confs.fillHeads();

        auto lnk = createLinkedCellParticles();
        lnk->add(confs.dataTensor(),confs.getTags(),{0,nBeads-1},{0,N-1});

        auto acc = std::make_shared<pimc::linkedCellAccelerationStructure>(lnk);

        confs.setAccelerationStructure(acc);


        return confs;
    }


    auto createGriddedPotentialAction()
    {
        auto G = std::make_shared<pimc::caoBernePropagator>(getTimeStep(),scatteringLength);

        auto meshKernel=std::make_shared<pimc::pairProductMeshKernel<pimc::caoBernePropagator> >(G,N,nBeads);

        meshKernel->setGeometry( createGeometry());

        auto SG = std::make_shared<pimc::actionTwoBodyMesh>();
        SG->setKernel(meshKernel);
        SG->setTimeStep(getTimeStep() );
        SG->setGeometry( createGeometry() );

        return SG;

    }


    
    auto createPotentialAction()
    {
        auto sTwoBodyCreator = std::make_shared<pimc::actionTwoBodyConstructor >() ;
        sTwoBodyCreator->setTimeStep(getTimeStep()  );
        sTwoBodyCreator->setGeometry( createGeometry());
        sTwoBodyCreator->setNBeads(nBeads);

        sTwoBodyCreator->setNMaxParticles(N);

        sTwoBodyCreator->registerGreenFunction<pimc::caoBernePropagatorTruncated>("caoBerne");

        
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

        
        return sG;
    }

    auto createAction( std::string kind)
    {
        std::shared_ptr<pimc::action> SG;
        if (kind== "exact" )
        {
            SG = createPotentialAction(); 
        }
        else if (kind == "gridded")
        {
            SG = createGriddedPotentialAction();
        }


        std::shared_ptr<pimc::action> sT= std::make_shared<pimc::kineticAction>(getTimeStep(), N , getNBeads()  ,  createGeometry() );

        auto S = std::make_shared<pimc::firstOrderAction> (sT,  SG);

        return S;
    }







/* 


    void run()
    {

     
        




        int seed=567;

        configurations.setRandom( lBox , randG );

        configurations.fillHeads();

        
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

    void SetDoCutOff( bool doCutOff)  {_doCutOff=doCutOff;}  */
    
    private:

    std::array<size_t,getDimensions()> nCells;
    std::array<Real,getDimensions()> lBox;
    Real beta;
    int nBeads;
    int N;
    Real scatteringLength;

    Real p;
    Real p2;
    randomGenerator_t randG;
    bool _doCutOff;
    int nBurnsIn;

};


/* TEST_F( caoBerneGridded , nCells3)
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

} */
/* TEST_F( caoBerneGridded , cutOff)
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

} */

/* TEST_F( caoBerneGridded , accuracy)
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
 */

TEST_F( caoBerneGridded , grid_allocation)
{
    int N=500,nBeads=10; 
    Real L = 1;


    SetN(N);
    SetNBeads(nBeads);
    SetBeta(1);
    SetBox({L,L,L});
    SetNCells( {6 ,6 ,6 } );
    SetScatteringLength(0.01);

    randomGenerator_t randG(567);

    auto configurations = createConfigurations( randG );

    configurations.fillHeads();
    auto link = createLinkedCellParticles();

    link->add(configurations.dataTensor(), {0,nBeads-1} , {0,N-1} );

    const auto & data = configurations.dataTensor();


    for(int t=0;t<=nBeads;t++)
    {
        for(int i=0;i<N;i++)
            for(int d=0;d<getDimensions();d++)
            {
                auto & cell = (*link)[t].getParticleCell(i);
                auto index= cell.getCellIndex();

                auto xd=pimc::pbc( data(i,d,t) , L , 1/L );

                // is the particle allocated to the right cell ?
                ASSERT_EQ( index[d] , (int)((xd+L/2 )/(L/getNCells()[d]) ) ) ;    

                // is the particle index stored in the cell ?
                const auto & particleIndices = cell.getParticleIndices();
                ASSERT_FALSE(   (std::find( particleIndices.begin(),particleIndices.end(), i) == particleIndices.end() )  );

                const auto & neighbours = cell.neighbours();

                ASSERT_EQ( neighbours.size() , 26 );




            
            }
    }





}


TEST_F( caoBerneGridded , evaluation)
{
    int N=500,nBeads=10; 
    Real L = 1;


    SetN(N);
    SetNBeads(nBeads);
    SetBeta(1);
    SetBox({L,L,L});
    SetNCells( {6 ,6 ,6 } );
    SetScatteringLength(0.01);
    randomGenerator_t randG(567);
    

    auto configurations = createConfigurations(randG);
    configurations.fillHeads();



    auto SC = createPotentialAction();
    auto SG = createGriddedPotentialAction();


    int iParticle=10;

    auto sEv= SC->evaluate(configurations,{0,nBeads-1}, iParticle);

    auto sEvG= SG->evaluate(configurations, {0,nBeads-1}, iParticle);





    for(int t=0;t<nBeads;t++)
    {
        auto n1=countDistances( configurations.dataTensor() , configurations.getAccelerationStructure().getLinkedCellList().cellLength(0),t,iParticle,{0,N-1}, createGeometry()  );

        auto kernel = static_cast<pimc::pairProductMeshKernel<pimc::caoBernePropagator>  & >( SG->getKernel() ); 
        auto n2 = kernel.getDistanceBuffer().getNDistances(t);
        ASSERT_EQ(n1 , n2);

    }

    ASSERT_NEAR(sEv,sEvG,TOL);

    sEv= SC->evaluate(configurations);
    sEvG= SG->evaluate(configurations);

    ASSERT_NEAR(sEv,sEvG,TOL);

    sEv= SC->evaluateTimeDerivative(configurations);
    sEvG= SG->evaluateTimeDerivative(configurations);



    ASSERT_NEAR(sEv,sEvG,TOL);

    Eigen::Tensor<Real,3> force( N , getDimensions() ,nBeads + 1);
    force.setConstant(0);
    Eigen::Tensor<Real,3> force2( N , getDimensions() ,nBeads + 1);
    force.setConstant(0);
    force2.setConstant(0);

    SC->addGradient(configurations,force);
    SG->addGradient(configurations,force2);

    for(int t=0;t<=nBeads;t++)
    {
        for(int i=0;i<N;i++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    ASSERT_NEAR( force(i,d,t),force2(i,d,t),TOL);
                }
            }
    }


    auto S1 = createAction("exact");
    auto S2= createAction("gridded");

    pimc::thermodynamicEnergyEstimator energyEstimator;
    pimc::virialEnergyEstimator virialEnergy(N, getNBeads());


    Real e1=energyEstimator(configurations,*S1);
    Real e2=energyEstimator(configurations,*S2);

    ASSERT_NEAR(e1,e2,TOL);
    
    Real eV1=virialEnergy(configurations,*S1);
    Real eV2=virialEnergy(configurations,*S2);

    ASSERT_NEAR(eV1,eV2,TOL);


    configurations.setHeadTail( iParticle, 8 , -1 );

    configurations.update( {0,nBeads-1}, {iParticle,iParticle} );
    sEv= SC->evaluate(configurations,{0,nBeads-1}, iParticle);
    sEvG= SG->evaluate(configurations, {0,nBeads-1}, iParticle);

    ASSERT_NEAR( sEv,sEvG,TOL);

    sEv= SC->evaluate(configurations);
    sEvG= SC->evaluate(configurations);

    ASSERT_NEAR(sEv,sEvG,TOL);

};


TEST_F( caoBerneGridded , update)
{
    int N=500,nBeads=10; 
    Real L = 1;


    SetN(N);
    SetNBeads(nBeads);
    SetBeta(1);
    SetBox({L,L,L});
    SetNCells( {5 ,5  ,5  } );

    SetScatteringLength(0.01);

    randomGenerator_t randG(567);
    auto configurations = createConfigurations(randG);
    configurations.fillHeads();

    auto SC = createPotentialAction();
    auto SG = createGriddedPotentialAction();

    int iParticle=10;

    auto sEv= SC->evaluate(configurations,{0,nBeads-1}, iParticle);

    auto sEvG= SG->evaluate(configurations, {0,nBeads-1}, iParticle);


    ASSERT_NEAR(sEv,sEvG,TOL);

    auto & data = configurations.dataTensor();

    int t0=0;
    int t1=nBeads - 1;
    std::normal_distribution<Real> gauss(0,0.2);

    for (int t=t0;t<=t1;t++)
    {
        for(int i=iParticle -2 ; i <= iParticle + 2;i++)
            for(int d=0;d<getDimensions();d++)
                {
                    data(i,d,t)+= L*gauss( randG );
                }   
    }

    configurations.fillHeads();

    sEv= SC->evaluate(configurations,{0,nBeads-1}, iParticle);

    configurations.update( {t0,t1},{iParticle-2,iParticle+2});
    sEvG= SG->evaluate(configurations, {0,nBeads-1}, iParticle);
    ASSERT_NEAR(sEv,sEvG,TOL);

};


TEST_F( caoBerneGridded , moves)
{
    int N=20,nBeads=20; 
    Real L = 4.940316131974357;
    
    SetN(N);
    SetNBeads(nBeads);
    SetBeta(2);
    SetBox({L,L,L});
    SetNCells( {3 ,3 ,3  } );
    Real a=0.0844782175465501;
    SetScatteringLength(a);
    randomGenerator_t randG1(567);
    randomGenerator_t randG2(567);

    auto configurations1 = createConfigurations(randG1);
    auto configurations2 = createConfigurations(randG2);

    auto SC = createAction("exact");
    auto SG = createAction("gridded");
    int l =0.6*getNBeads();

    pimc::levyMove levy1(l,0);
    pimc::levyMove levy2(l,0);
    

    pimc::thermodynamicEnergyEstimator energyEstimator;
    pimc::virialEnergyEstimator virialEnergy(N, getNBeads());

    int nTrials=1000;
    Real sum1=0,e1=0,eV1=0;
    Real sum2=0,e2=0,eV2=0 ;
    for(int iTrial=0;iTrial < nTrials ; iTrial ++ )
    {
        levy1.attemptMove(configurations1,*SC,randG1);
        sum1= SC->getPotentialAction().evaluate(configurations1);
        e1=energyEstimator(configurations1,*SC);
        eV1=virialEnergy(configurations1,*SC);

        levy2.attemptMove(configurations2,*SG,randG2);
        configurations2.update({0,nBeads-1},{0,N-1});
        sum2= SG->getPotentialAction().evaluate(configurations2);
        e2=energyEstimator(configurations2,*SG);
        eV2=virialEnergy(configurations2,*SG);

        ASSERT_NEAR( sum1,sum2,TOL);

        ASSERT_NEAR( e1,e2,TOL);
        ASSERT_NEAR( eV1,eV2,TOL); 
        
    }



}