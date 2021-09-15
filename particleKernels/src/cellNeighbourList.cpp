#include "cellNeighbourList.h"
#include <array>
#include <cmath>
#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <memory.h>


using namespace particleKernels;

template<int d>
inline int indexFortranStorage(const int* im, const int* sizes)
{
    return sizes[0]*(  sizes[1]*im[2] + im[1] )  + im[0] ;
}


template<int d>
void indicesFromIndexFortranStorage(int* im, int cellIndex,const int* sizes)
{
    im[0]=cellIndex%sizes[0];
    cellIndex-=im[0];
    cellIndex/=sizes[0];


    im[1]=cellIndex%sizes[1];
    cellIndex-=im[1];
    cellIndex/=sizes[1];


    im[2]=cellIndex;



}


template<class allocator_t, class T>
T* reallocate(allocator_t alloc, T* data , size_t oldSize, size_t newSize)
{
    T* newData=alloc.allocate(newSize);
    memcpy(newData, data, oldSize * sizeof(T));

    alloc.deallocate(data,oldSize);

    return newData;

}



int gridIndexPerDimensionRestrictedToBox( double x , double lBox, int nCell)
{
    return std::floor( (  x/lBox + 0.5 ) *nCell);
};

inline int gridIndexPerDimensionUnRestrictedToBox( double x , double lBox, int nCell)
{
    int iBox = std::round( x/lBox );
    
    return gridIndexPerDimensionRestrictedToBox(x-iBox*lBox,lBox,nCell);

};

int gridIndexPerDimension( double x , double lBox, int nCell)
{
    return gridIndexPerDimensionRestrictedToBox(x,lBox,nCell);

};

 


template<int dimensions, class allocator_t=std::allocator<int>  >
void listFortranStorage(
    const Real * positions, 
    const Real* lBox, const int * nCells, 
    int** cells, int * nParticlesPerCell,int * nParticlesPerCellAllocated,
    int * cellIndexPerParticle,   int * subIndexPerParticle,
     allocator_t & alloc,
     int startParticle,int endParticle, int N)
{

    constexpr int buffer = 10;


    std::array<int ,dimensions> indexMultiple;



    for(int i=startParticle;i<=endParticle;i++)
    {
        for(int d=0;d<dimensions;d++)
        {

            indexMultiple[d]=gridIndexPerDimension(  positions[i + d *N],lBox[d] ,nCells[d] );


            assert(indexMultiple[d] >=0 & indexMultiple[d] < nCells[d]  );
        }


        int index=indexFortranStorage<dimensions>(indexMultiple.data() ,nCells  );

        cellIndexPerParticle[i]=index;


        // reallocate if required , should not happen very often except at the beginning of the simulation
        if ( nParticlesPerCellAllocated[index] <= nParticlesPerCell[index] )
        {
            
            cells[index]=reallocate(alloc, cells[index],nParticlesPerCell[index]
            , nParticlesPerCell[index] + buffer );

            nParticlesPerCellAllocated[index]=nParticlesPerCell[index] + buffer;



        }

        // record the position of the particle
        int * particlesInCell=cells[index];
        particlesInCell[ nParticlesPerCell[index] ]=i;


        subIndexPerParticle[i]=nParticlesPerCell[index];
        nParticlesPerCell[index]+=1;


    }

}

template<int dimensions, class allocator_t=std::allocator<int>  >
void updateListFortranStorage(
    const Real * positions, 
    const Real* lBox, const int * nCells, 
    int** cells, int * nParticlesPerCell,int * nParticlesPerCellAllocated,
    int * cellIndexPerParticle,   int * subIndexPerParticle,
     allocator_t & alloc,
     int startParticle,int endParticle, int N)
{

    constexpr int buffer = 10;

    std::array<int ,dimensions> indexMultiple;


    for(int i=startParticle;i<=endParticle;i++)
    {
        for(int d=0;d<dimensions;d++)
        {

            indexMultiple[d]=gridIndexPerDimension(  positions[i + d *N],lBox[d] ,nCells[d] );


            assert(indexMultiple[d] >=0 & indexMultiple[d] < nCells[d]  );
        }


        int index=indexFortranStorage<dimensions>(indexMultiple.data() ,nCells  );

        if ( cellIndexPerParticle[i]!=index )
        {
            // removes the particle from the old grid cell
            
            int iCellOld=cellIndexPerParticle[i];
            int iiCellOld=subIndexPerParticle[i];

            std::swap( (cells[iCellOld])[iiCellOld],(cells[iCellOld])[
                nParticlesPerCell[iiCellOld] - 1
            ] );

            nParticlesPerCell[iiCellOld]-=1;

            // assign the particle to the new grid cell
            cellIndexPerParticle[i]=index;

            

            // reallocate if required , should not happen very often except at the beginning of the simulation
            if ( nParticlesPerCellAllocated[index] <= nParticlesPerCell[index] )
            {
                
                cells[index]=reallocate(alloc, cells[index],nParticlesPerCell[index]
                , nParticlesPerCell[index] + buffer );

                nParticlesPerCellAllocated[index]=nParticlesPerCell[index] + buffer;
            }


            
            int * particlesInCell=cells[index];
            particlesInCell[ nParticlesPerCell[index] ]=i;


            subIndexPerParticle[i]=nParticlesPerCell[index];
            nParticlesPerCell[index]+=1;

            }

    }

}





void simpleCellNeighbourList::list(const double * positions, size_t N, int iStart, int iEnd)
{

    if (N > maxNParticles)
    {

        cellIndexPerParticle=reallocate(intAllocator,cellIndexPerParticle,maxNParticles,N);

        subIndexPerParticle=reallocate(intAllocator,subIndexPerParticle,maxNParticles,N);

    }
    if ( dimensions==3)
    {
        listFortranStorage<3>(positions,lBox.data(),nCells.data(),cells.data() , nParticlesPerCell.data(),nParticlesPerCellAllocated.data() ,cellIndexPerParticle,subIndexPerParticle,
        intAllocator,iStart,iEnd,N);
    }
    else
    {
        throw std::invalid_argument("List with dimensions !=3 not supported");


    }

}

void simpleCellNeighbourList::updateList(const double * positions, size_t N, int iStart, int iEnd)
{

    if ( dimensions==3)
    {
        updateListFortranStorage<3>(positions,lBox.data(),nCells.data(),cells.data() , nParticlesPerCell.data(),nParticlesPerCellAllocated.data() ,cellIndexPerParticle,subIndexPerParticle,
        intAllocator,iStart,iEnd,N);
    }

}




simpleCellNeighbourList::simpleCellNeighbourList(std::vector<int> nCells_, std::vector<Real> lBox_, neighbourListPolicy policy_,  int nParticleStorageInitial ) :
dimensions(nCells_.size() ),
policy(policy_),
nCells(nCells_),
lBox(lBox_),
periodic(true),
maxNParticles(10000)
{
    nCellsTotal=1;
    lBoxCell.resize(3);

    cellIndexPerParticle=intAllocator.allocate(maxNParticles);
    subIndexPerParticle=intAllocator.allocate(maxNParticles);


    for(int d=0;d<dimensions;d++)
    {
        nCellsTotal*=nCells[d];
        lBoxCell[d]=lBox[d]/nCells[d];
    }

    cells.resize(nCellsTotal,nullptr);
    indexCellNeighbours.resize(nCellsTotal,nullptr);
    displacementCellNeighbours.resize(nCellsTotal,nullptr);

    nNeighboursPerCell.resize(nCellsTotal,0);

    nParticlesPerCell.resize(cells.size() ,0);
    nParticlesPerCellAllocated.resize(cells.size(),0);

    nNeighboursPerCellMax=std::pow(3,dimensions) ;

    // allocate acceleration data structure
    for (int i=0; i<cells.size();i++)
    {
        cells[i]=intAllocator.allocate(nParticleStorageInitial);
        nParticlesPerCell[i]=0;
        nParticlesPerCellAllocated[i]=nParticleStorageInitial;    
        indexCellNeighbours[i]=intAllocator.allocate(  nNeighboursPerCellMax  );
        indexCellNeighbours[i]=intAllocator.allocate(  nNeighboursPerCellMax  );
        displacementCellNeighbours[i]=realAllocator.allocate(  nNeighboursPerCellMax * dimensions  );

        for (int ii=0;ii<nNeighboursPerCellMax*dimensions;ii++)
        {
            displacementCellNeighbours[i][ii]=0;
        }


    }

    // store indices of neighbouring clusters

        if (dimensions == 3 )
        {
            for (int i=0;i<nCells[0];i++)
            {
                for (int j=0;j<nCells[1];j++)
                {
                    for (int k=0;k<nCells[2];k++)
                    {
                        std::array<int,3> mIndex{i,j,k};

                        int cellIndex =indexFortranStorage<3>(mIndex.data(),nCells.data());

                        assert( cellIndex < nCellsTotal);


                        // loop over neighbouring indices
                        int iNeighbour=0;
                        for (int ii=i-1;ii<=i+1;ii++)
                        {
                            for (int jj=j-1;jj<=j+1;jj++)
                                for (int kk=k-1;kk<=k+1;kk++)
                                {
                                                                                                               int iii=ii;
                                    int jjj=jj;
                                    int kkk= kk;


                                    if ( periodic)
                                    {
                                        // low face
                                        if (ii < 0 )
                                        {
                                            iii+=nCells[0];
                                            displacementCellNeighbours[cellIndex][iNeighbour]=-lBox[0];
                                        }

                                        if (jj < 0 )
                                        {
                                            jjj+=nCells[1];
                                            displacementCellNeighbours[cellIndex][iNeighbour + nNeighboursPerCellMax]=-lBox[1];
                                        }

                                        if (kk < 0 )
                                        {
                                            kkk+=nCells[2];
                                            displacementCellNeighbours[cellIndex][iNeighbour + 2*nNeighboursPerCellMax]=-lBox[2];
                                        }

                                        // high face

                                        if (ii >=nCells[0] )
                                        {
                                            iii-=nCells[0];
                                            displacementCellNeighbours[cellIndex][iNeighbour]=+lBox[0];
                                        }

                                        if (jj >= nCells[1] )
                                        {
                                            jjj-=nCells[1];
                                            displacementCellNeighbours[cellIndex][iNeighbour + nNeighboursPerCellMax]=lBox[1];
                                        }

                                        if (kk >= nCells[2] )
                                        {
                                            kkk-=nCells[2];
                                            displacementCellNeighbours[cellIndex][iNeighbour + 2*nNeighboursPerCellMax]=lBox[2];
                                        }





                                    }


                                   std::array<int,3> mIndexCell={iii,jjj,kkk};
                                    int cellIndexNeighbour =indexFortranStorage<3>(mIndexCell.data(),nCells.data());

                                    

                                        
                                        





                                            if (
                                                (iii>=0 and iii<nCells[0] ) and(jjj>=0 and jjj<nCells[1] ) and 
                                                (kkk>=0 and kkk<nCells[2] ) 
                                                ) 
                                            {

                                                
                                                assert(cellIndexNeighbour >= 0 and cellIndexNeighbour<nCellsTotal  );
                                    
                                                // store the neighbours
                                                indexCellNeighbours[cellIndex][iNeighbour]=cellIndexNeighbour;       iNeighbour++;

                                                nNeighboursPerCell[cellIndex]=iNeighbour;
                                            }

                                        

                                                              

                                    

                                 

                                }
                            

                        }







                    }
                }
            }

        }
        else
        {
            throw std::invalid_argument("Dimensions!=3 not supported");
        }

}

simpleCellNeighbourList::~simpleCellNeighbourList()
{
    for (int i=0; i<cells.size();i++ )
    {
        intAllocator.deallocate(cells[i],nParticlesPerCellAllocated[i]);
        intAllocator.deallocate(indexCellNeighbours[i],nNeighboursPerCellMax);

        

    }
    intAllocator.deallocate(cellIndexPerParticle,maxNParticles);


}


std::string simpleCellNeighbourList::reportCellNeighbourTopology() const
{
    std::stringstream reportStream;
    std::string delim = "\t";

    reportStream << "i" << delim << "j" << delim << "k" << delim << "cellIndex" ; 
    for (int iNeigh=0;iNeigh<nNeighboursPerCellMax; iNeigh++)
    {
        reportStream << delim << "neighBour" << iNeigh ;
    }
    reportStream << std::endl;

    reportStream << "-------------------------------" << std::endl;
    if ( dimensions == 3 )
    {
        for (int i=0;i<nCells[0];i++)
            for (int j=0;j<nCells[1];j++)
                for (int k=0;k<nCells[2];k++)

                {
                    std::array<int,3> mIndex= {i,j,k};
                    int cellIndex=indexFortranStorage<3>( mIndex.data(),nCells.data()) ;


                    reportStream << i << delim << j << delim << k << delim << cellIndex;

                    for (int iNeigh=0;iNeigh<nNeighboursPerCell[cellIndex];iNeigh++)
                    {
                        reportStream << delim << indexCellNeighbours[cellIndex][iNeigh]; 
                    }
                    for (int iNeigh=nNeighboursPerCell[cellIndex];iNeigh<nNeighboursPerCellMax; iNeigh++)
                    {
                        reportStream << delim << "-";
                    }


                    reportStream<< std::endl;




                }

    }
    else
    {
        throw std::invalid_argument("Dimensions != 3 not yet supported");
    }


    return reportStream.str();

}

void simpleCellNeighbourList::sanityCheck()
{
    assert(nCellsTotal == cells.size() );

    int nMaxNeighbours=std::pow(3,dimensions);

    for (int i=0;i<cells.size();i++)
    {

        std::array<int,3> mIndexCell;
         indicesFromIndexFortranStorage<3>(mIndexCell.data(), i, nCells.data() );

        assert(nNeighboursPerCell[i]<=nNeighboursPerCellMax);

        std::array<int,3> nWrapped{0,0,0};
        std::array<int,3> nWrappedExpected{0,0,0};

        if (periodic)
        {
            for (int d=0;d<dimensions;d++)
            {
                if ( mIndexCell[d]==0  )
                {
                    nWrappedExpected[d]+=std::pow(3,dimensions-1);
                }

                if ( mIndexCell[d]==nCells[d] -1  )
                {
                    nWrappedExpected[d]+=std::pow(3,dimensions-1);
                }



            }

            assert(nNeighboursPerCell[i]==std::pow(3,dimensions));
        }

        for (int ii=0;ii<nNeighboursPerCell[i];i++)
        {
            
            int cellNeighbour=indexCellNeighbours[i][ii];


            if ( dimensions == 3 and periodic )
            {
                std::array<int,3> mIndexCellN;

                 indicesFromIndexFortranStorage<3>(mIndexCellN.data(), cellNeighbour, nCells.data() );

                 for (int d=0;d<dimensions;d++)
                 {
                     if ( ( mIndexCellN[d] == (nCells[d] - 1) ) and
                        mIndexCell[d] == 0
                       )
                     {

                         assert(displacementCellNeighbours[i][ii + d*nMaxNeighbours] == -lBox[d]);

                         nWrapped[d]+=1;
                     }

                     if ( ( mIndexCell[d] == nCells[d] - 1 ) and
                        mIndexCellN[d] == 0
                       )
                     {
                         assert(displacementCellNeighbours[i][ii + d*nMaxNeighbours] == lBox[d]);
                         nWrapped[d]+=1;
                     }



                 }

            }

           
           

            

            // check if the current cell is listed in the neighbour cells list of cellNeighbour

            bool foundCurrentCell=false;

            for (int jj=0;jj<nNeighboursPerCell[cellNeighbour];jj++)
            {
                assert(
                    ( indexCellNeighbours[cellNeighbour][jj] >= 0 )
                and 
                ( indexCellNeighbours[cellNeighbour][jj] >= 0 )


                );
                if (indexCellNeighbours[cellNeighbour][jj]==i)
                foundCurrentCell=true;

            }

            assert(foundCurrentCell==true);

        }



         for (int d=0;d<dimensions;d++)
            {
                
                assert(nWrappedExpected[d]==nWrapped[d]);
            }


    }


}


  std::string simpleCellNeighbourList::reportParticleDistribution(const Real * positions, int iStart, int iEnd, int N) const
    {
        const std::string delim="\t";
        std::stringstream ss;
        std::vector<int> mIndex{0,0,0};


        ss << "particle" << delim << "cell";

        for ( int d=0;d<dimensions;d++)
        {
            ss << delim << "X[" << d << "]";
        }

        ss << std::endl;

        for (int i=iStart;i<=iEnd;i++)
        {
            auto iCell = cellIndexPerParticle[i];
            auto iiCell= subIndexPerParticle[i];

            assert(    (cells[iCell])[iiCell] == i  );
            
            
            if (dimensions==3)
            {
                indicesFromIndexFortranStorage<3>(mIndex.data(), iCell,
                
                nCells.data() );

            }
            else
            {
                throw std::invalid_argument("D!=3 not supported");
            }            


            ss << i << delim << iCell ;

             for (int d=0; d<dimensions;d++)
            {
                ss << delim <<  mIndex[d];

           

                 int indexDCheck=std::floor( (positions[i + d*N] / lBox[d] + 0.5 )*nCells[d] );

                 assert(mIndex[d] == indexDCheck  );
            }


            for (int d=0; d<dimensions;d++)
            {
                ss << delim <<  positions[i + d*N];
            }
            
            ss << std::endl ;
        }
        
        return ss.str();


    }

