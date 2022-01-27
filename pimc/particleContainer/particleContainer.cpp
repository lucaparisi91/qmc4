

#include "particleContainer.h"
#include <array>
#include <cmath>
#include <cassert>
#include <exception>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <memory.h>

namespace pimc
{

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

auto wrap( int index, int N)
{
    if (index < 0) {return index + N;};
    if (index >= N) {return index - N ;}

    return index;
};

auto nDisplacement( int index , int N )
{
    if (index<0) return +1;
    if (index>=N) return - 1;

    return 0;
    
}

simpleCellNeighbourList::simpleCellNeighbourList(std::array<size_t,getDimensions()> nCells, std::array<Real,getDimensions()> lBox   ) :
_nCells(nCells),
_lBox(lBox),
nParticles(0)
{
    int nCellsTotal=1;
    for(int d=0;d<getDimensions();d++)
    {
        nCellsTotal*=nCells[d];
        _left[d]=-lBox[d]/2;
        _right[d]=lBox[d]/2;
        _delta[d]=_lBox[d]/_nCells[d];
    }

    cells.resize(nCellsTotal);


    for ( size_t i=0;i<nCellsTotal;i++)
    {
        cells[i]=(new cell{});
    }


    for (int i=0;i<_nCells[0];i++)
        for (int j=0;j<_nCells[1];j++)
            for (int k=0;k<_nCells[2];k++)
                    {            
                        auto & neighbours = getCell(i,j,k).neighbours();
                        auto & displacements = getCell(i,j,k).displacements();

                        getCell(i,j,k).setCellIndex({i,j,k});


                        int iNeigh=0;
                        for (int ii=i-1;ii<=i+1;ii++)
                            for (int jj=j-1;jj<=j+1;jj++)
                                for (int kk=k-1;kk<=k+1;kk++)
                                {

                                    if ( (ii==i) and (jj==j) and (kk==k) )
                                    { 
                                        continue;
                                    }

                                    int iii=wrap(ii,_nCells[0]);
                                    int jjj=wrap(jj,_nCells[1]);
                                    int kkk=wrap(kk,_nCells[2]);

                                    neighbours[iNeigh]= & (getCell(iii,jjj,kkk) );
                                    
                                    displacements(iNeigh,0)=nDisplacement(ii,_nCells[0])*_lBox[0];
                                    displacements(iNeigh,1)=nDisplacement(jj,_nCells[1])*_lBox[1];
                                    displacements(iNeigh,2)=nDisplacement(kk,_nCells[2])*_lBox[2];

                                    iNeigh++;
                                }      
                     }
}


simpleCellNeighbourList::~simpleCellNeighbourList()
{
   for ( int i=0;i<cells.size();i++)
   {
       delete cells[i];
   }
}

void simpleCellNeighbourList::setCapacity(size_t N)
{

    cellIndexPerParticle.resize( N );
    subIndexPerParticle.resize( N);


    for ( auto cell : cells)
    {
        cell->setCapacity(N);
    }


}


cell::cell( size_t nMaxParticles) : _nMaxParticles(nMaxParticles),
    _positions( nMaxParticles,getDimensions( )),
    _displacements( std::pow(3,getDimensions() ) - 1 ,getDimensions()   ),
    _nParticles(0),
    _index{0,0,0} {
        int nNeighbours=std::pow(3,getDimensions()) -1;

        _displacements.setConstant(0);
        _neighbours.resize( nNeighbours,nullptr);
    };

void cell::setCapacity(size_t N)
{
    _positions.resize( N, getDimensions());
    _particleIndex.resize(N);
    _nMaxParticles=N;
}


linkedCellParticles::linkedCellParticles( std::array<size_t,getDimensions()> nCells, std::array<Real,getDimensions() > lBox  ) :
_nCells(nCells),
_lBox(lBox)
 {}

void linkedCellParticles::setCapacity(size_t N,size_t M)
{
    auto oldSize=particles.size();
    particles.resize( M + 1);
    for(int i=oldSize;i<M + 1;i++)
    {
        particles[i]=std::make_shared<linkedCell_t>(_nCells,_lBox);
    };

    for(int i=0;i<M + 1;i++)
    {
        particles[i]->setCapacity(N);
    }

}

void linkedCellParticles::add( const Eigen::Tensor<Real,3> & data, const range_t &  timeRange,const  range_t & particleRange)
{

    for (int t=timeRange[0]; t<=timeRange[1]+1;t++)
    {
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {
            particles[t]->add(i, data(i,0,t),data(i,1,t),data(i,2,t));
        }
    }
}

void linkedCellParticles::remove(  const range_t &  timeRange,const  range_t & particleRange)
{

    for (int t=timeRange[0]; t<=timeRange[1];t++)
    {
        for(int i=particleRange[0];i<=particleRange[1];i++)
        {
            particles[t]->remove(i);
        }
    }
}






}


