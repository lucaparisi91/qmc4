#include <vector>
#include <array>
#include <cmath>
#include <string>
#include <string>
#include "traits.h"
#include "pimcConfigurations.h"

namespace acceleration
{

using Real = double;

class neighbourListPolicy
{};


class particlesRestrictedInBoxPolicy :public neighbourListPolicy
{};



class simpleCellNeighbourList
{

    public:
    /*
    Naive data structure acceleration for short ranged ranged interactions.
    Subdivide the box in nCells[d] cell per dimension d. 
    Stores acceleration data structure as an array of pointers to contigous block of data. The blocks are reallocated on the run if required
    */

   using realAllocator_t = std::allocator<Real> ;
   using intAllocator_t = std::allocator<int> ;


    simpleCellNeighbourList(std::array<int,DIMENSIONS >  nCells, std::array<Real,DIMENSIONS> lBox, neighbourListPolicy = particlesRestrictedInBoxPolicy() ,  int nParticleStorageInitial = 10   );


    void add(const double * positions, size_t N, int iStart, int iEnd); // inserts particle in the cell neighbour list

    void update(const double * positions, size_t N, int iStart, int iEnd); // update the cell index of the particle

    void remove(  int iStart, int iEnd ); // remove the particle from the partition structure

    ~simpleCellNeighbourList();

    void sanityCheck();


    std::string reportCellNeighbourTopology() const;

    std::string reportParticleDistribution(const Real * positions, int iStart, int iEnd , int N) const;


    inline int cellIndex(int iParticle) const {return cellIndexPerParticle[iParticle];}

    int nCellNeighbours(int iCell) const  { return nNeighboursPerCell[iCell];}

    auto neighbourCell( int iCell,  int iiCell) const {
        return indexCellNeighbours[iCell][iiCell];
    }



    int  nParticles( int iCell) const {
        return nParticlesPerCell[iCell];
    }

    auto particle(int iCell, int iParticle ) const
    {
        return cells[iCell][iParticle];
    }


    auto nParticles() const { return nTotParticles; }

    private:
 
    inline int index(int iX, int iY, int iZ) const {return iZ* nCells[0]*( nCells[1] + iY )  + iX  ;} // transform a 3d index in a 1d intex
    inline int index(int iX, int iY) const {return  nCells[0]* iY  + iX  ;} // transform a 2d index in a 1d intex
    inline int index(int iX) const {return  iX  ;} // transform a 1d index in a 1d intex

    inline int index(const std::array<int,3> & im) const {return index(im[0],im[1],im[2]);}

    int dimensions;
    realAllocator_t realAllocator;
    intAllocator_t intAllocator;

    // subdivision properties
    std::array<Real,DIMENSIONS> lBoxCell;
    std::array<Real,DIMENSIONS> lBox;
    std::array<int,DIMENSIONS> nCells;

    // data structure to contain the decomposition data in cells
    std::vector<int*> cells; // one entry per cell
    std::vector<int> nParticlesPerCell; // number of particles in a certain cell
    std::vector<int> nParticlesPerCellAllocated; // number of particles allocated in a certain cell
    
    /*  
    Store grid information associated with ich particles
    len(cellIndexPerParticle)>= Nmax particles
    cells[ cellIndexPerParticle[i] ][ subIndexPerParticle[i]   ] = i   */

    int * cellIndexPerParticle;
    int * subIndexPerParticle;

    int maxNParticles;

    std::vector<int*> indexCellNeighbours; // number of neighBours in the cell
    std::vector<int> nNeighboursPerCell;
    
    std::vector<Real*> displacementCellNeighbours; // contains the displacement to add when cells wrap around the boundary

    int nNeighboursPerCellMax;

    bool periodic;

    neighbourListPolicy policy;
    int nCellsTotal;
    int nParticlesTotal;
    int nTotParticles;

};

}


namespace pimc
{
    struct cellNeighbourList
    {
        using data_t = Eigen::Tensor<Real,3>;

        void setGeometry(const geometry_t & geo);
        void setNBeads(int nBeads_){_nBeads=nBeads_;};
        void setNCells(const std::array<int,DIMENSIONS > & nCells_);

        void build();

        void add(const data_t & configurations, const std::array< int, 2  > & timeRange , const std::array< int, 2  > & particleRange );

        void remove(const std::array< int, 2 > & timeRange , const std::array< int, 2  > & particleRange);

        void update(const data_t  & configurations, const std::array< int, 2  > & timeRange , const std::array< int, 2  > & particleRange);

        inline auto & cellNeighBourList( int t) {
            return *(_timeSliceLists[t]);
        }

        inline const auto & cellNeighBourList( int t) const {
            return *(_timeSliceLists[t]);
        }


        private:

        std::vector< std::unique_ptr<acceleration::simpleCellNeighbourList>  > _timeSliceLists;
        std::array<Real,DIMENSIONS> _lBox;
        std::array<int,DIMENSIONS> _nCells;
        int _nBeads;
};




};