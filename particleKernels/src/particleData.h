#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H



#include "lattice.h"
#include "cellData.h"
#include "cassert"

namespace particleKernels
{

template<class T>
class particleDataBase
{
    public:

    using cell_t = T;
    using index_t = int;

    particleDataBase(std::shared_ptr<lattice> _lattice);

    virtual void list(const Real * positions, int iStart, int iEnd, int N)=0;



    //virtual void updateList(const Real * positions, int iStart,int iEnd, int N)=0;

    const auto & getLattice() const {return *(_lattice);}

    int cellIndex(int iParticle) const  {return _cellIndexPerParticle[iParticle];}

    int subCellIndex(int iParticle) const {
        
        return _subCellIndexPerParticle[iParticle];        
     }

    
    void updateList(const double * positions, int iStart, int iEnd, int N);
    auto & operator[](size_t i) noexcept {return *_cellData[i];}
    const auto & operator[](size_t i) const noexcept {return *_cellData[i];}

   

    void resizeParticleData(size_t nMax)
    {
        _cellIndexPerParticle.resize(nMax,-1);
        _subCellIndexPerParticle.resize(nMax,-1);
    }


    void recordParticleCellPosition(index_t iParticle,index_t iCell , index_t iSubCell)
    {
        _cellIndexPerParticle[iParticle]=iCell;
        _subCellIndexPerParticle[iParticle]=iSubCell;
    }

    void checkPBCGhostCells();

    private:


    std::shared_ptr<lattice> _lattice;

    /* Per particle data */
    std::vector<int>  _cellIndexPerParticle;
    std::vector<int> _subCellIndexPerParticle;

    /* Cell data  */
    std::vector<std::shared_ptr<cell_t> > _cellData;
};


class particleDataIndex : public particleDataBase<indexCellData>
{
    public:

    particleDataIndex(std::shared_ptr<lattice> lattice) : particleDataBase<indexCellData>(lattice) {}

    void list(const double * positions, int iStart, int iEnd, int N);


    private:
};


class particleData3D : public particleDataBase<particleCellData3d>
{
    public:

    particleData3D(std::shared_ptr<lattice> lattice) : particleDataBase<particleCellData3d>(lattice) {}


    void list(const double * positions, int iStart, int iEnd, int N);

   

    private:

};

}

#endif