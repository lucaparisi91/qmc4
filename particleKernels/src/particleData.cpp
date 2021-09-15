#include "particleData.h"
#include <iostream>
#include "utils.h"

using namespace particleKernels;


template<class T>
particleDataBase<T>::particleDataBase(std::shared_ptr<lattice> lattice_) :
_lattice(lattice_)
{
    const auto & currentLattice = getLattice();

    const int latticeSize = currentLattice.extendedSize();

    for (int i=0;i<currentLattice.extendedSize();i++)
    {
        _cellData.emplace_back( std::make_shared<cell_t>()  );
    };


    /* make periodic ghost cells point to image cells */

    if ( currentLattice.dimensions() == 3 )
    {
        for (int i=0;i<(currentLattice.extendedShape())[0];i++)
            for (int j=0;j<currentLattice.extendedShape()[1];j++)
                for (int k=0;k<currentLattice.extendedShape()[2];k++)
                {
                    std::array<index_t,3> index3_2{i,j,k};

                    for (int d=0;d<currentLattice.dimensions();d++)
                    {
                        if (index3_2[d]<currentLattice.lowIndex(d))
                        {
                            index3_2[d]+= currentLattice.shape()[d];
                        }
                         if (index3_2[d]>currentLattice.highIndex(d))
                        {
                            index3_2[d]-= currentLattice.shape()[d];
                        }   
                    }
                
                    auto currentIndex=currentLattice.index(i,j,k);
                    auto imageIndex=currentLattice.index(index3_2);

                    if (currentIndex != imageIndex)
                    {
                        _cellData[currentIndex]=_cellData[imageIndex];
                    }
                    
                }
    }
}


void particleDataIndex::list(const Real * positions,int iStart,int iEnd,int N)
{
    const auto & currentLattice=getLattice();
    const auto & lBox = currentLattice.lengthBox();
    const auto & lBoxInverse = currentLattice.inverseLengthBox();
    const auto & left = currentLattice.lowerEdge();

    resizeParticleData(N);


    if (currentLattice.dimensions()== 3 )
        {
        for (int i=0;i<N;i++) 
        {
            
            size_t index = currentLattice.index(
                positions[i],
                positions[i+N],
                positions[i+2*N]
            );

            assert(currentLattice.extendedSize() > index);
            assert(index > currentLattice.nGhosts()[0]);



            auto ii = (*this)[index].add(i);

            recordParticleCellPosition(i, index, ii);
        }
    }

};



void particleData3D::list(const Real * positions,int iStart,int iEnd,int N)
{
    const auto & currentLattice=getLattice();
    const auto & lBox = currentLattice.lengthBox();
    const auto & lBoxInverse = currentLattice.inverseLengthBox();
    const auto & left = currentLattice.lowerEdge();

    resizeParticleData(N);


    if (currentLattice.dimensions()== 3 )
        {
        for (int i=0;i<N;i++) 
        {            
            size_t index = currentLattice.index(
                positions[i],
                positions[i+N],
                positions[i+2*N]
            );

            assert(currentLattice.extendedSize() > index);
            assert(index > currentLattice.nGhosts()[0]);

            auto ii = (*this)[index].add(i,positions[i],positions[i+N],positions[i+2*N]);

            assert ( 
                (*this)[index].getParticlePosition(ii,0) == positions[i]
            );

            recordParticleCellPosition(i, index, ii);
        }
    }

};



template<class T>
void particleDataBase<T>::updateList(const Real * positions,int iStart,int iEnd,int N)
{
    const auto & currentLattice=getLattice();
    const auto & lBox = currentLattice.lengthBox();
    const auto & lBoxInverse = currentLattice.inverseLengthBox();
    const auto & left = currentLattice.lowerEdge();


    if (currentLattice.dimensions()== 3 )
        {
        for (int i=iStart;i<=iEnd;i++) 
        {
            auto newIndex = currentLattice.index(
                positions[i],
                positions[i+N],
                positions[i+2*N]
            );

            auto oldIndex = cellIndex(i);
            auto oldSubIndex = subCellIndex(i);

            
            (*this)[oldIndex].updatePosition(oldSubIndex,positions[i],positions[i+N],positions[i+2*N]);


            if (oldIndex != newIndex)
            {
                
                auto ii = (*this)[oldIndex].moveTo( oldSubIndex, (*this)[newIndex]  );
                auto iSwappedParticle = (*this)[oldIndex].getParticleIndex(oldSubIndex); 
                recordParticleCellPosition(iSwappedParticle, oldIndex,oldSubIndex);
                recordParticleCellPosition(i, newIndex, ii);

            }

            
            
        }
        }
};




template<class T>
void particleDataBase<T>::checkPBCGhostCells()
{
    const auto & currentLattice = getLattice();

    if ( currentLattice.dimensions() == 3 )
    {
        for (int i=currentLattice.lowIndex(0);i<=currentLattice.highIndex(0);i++)
            for (int j=currentLattice.lowIndex(1);j<=currentLattice.highIndex(1);j++)
                for (int k=currentLattice.lowIndex(2);k<=currentLattice.highIndex(2);k++)
                {
                    auto iCell = currentLattice.index(i,j,k);
                    std::array<index_t,3> index3_1=currentLattice.index3d(iCell);
                    assert(index3_1[0]==i);
                    assert(index3_1[1]==j);
                    assert(index3_1[2]==k);


                    for (int iOffset=0;iOffset<currentLattice.nCellsNeighbourhood() ; iOffset++)
                    {
                        auto iCell2 = currentLattice.getNeighbour(iCell,iOffset);


                        auto index3_2 = currentLattice.index3d(iCell2);
                        auto index3_2Old = currentLattice.index3d(iCell2);


                        for (int d=0;d<currentLattice.dimensions();d++)
                    {

                        if (index3_2[d]<currentLattice.lowIndex(d))
                        {
                            index3_2[d]+= currentLattice.shape()[d];
                        }
                         if (index3_2[d]>currentLattice.highIndex(d))
                        {
                            index3_2[d]-= currentLattice.shape()[d];
                        }   
                    }
                
                    auto imageIndex=currentLattice.index(index3_2);


                    if (iCell2 != imageIndex)
                    {
                       // std::cout << _cellData[iCell2]->nParticles() << " " << _cellData[imageIndex]->nParticles() << " " << index3_2Old[0] << " " << index3_2Old[1] << " " << index3_2Old[2] << " " << index3_2[0] << " " << index3_2[1] << " " << index3_2[2] << std::endl;
                        
                        
                        assert( _cellData[iCell2]->nParticles()==_cellData[imageIndex]->nParticles() );
                    }



                    }
                    
                    
                    
                }
    }
}


//template class particleDataBase<indexCellData>;
template class particleDataBase<particleCellData3d>;







