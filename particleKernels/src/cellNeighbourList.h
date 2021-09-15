

#include <vector>
#include <array>
#include <cmath>

namespace particleKernels
{

using Real = double;

class neighbourListPolicy
{};

class particlesRestrictedInBoxPolicy :public neighbourListPolicy
{};

#include "listKernels.hpp"


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


    simpleCellNeighbourList(std::vector<int> nCells, std::vector<Real> lBox, neighbourListPolicy = particlesRestrictedInBoxPolicy() ,  int nParticleStorageInitial = 10   );

    void list(const double * positions, size_t N, int iStart, int iEnd); // inserts particle in the cell neighbour list

    void updateList(const double * positions, size_t N, int iStart, int iEnd); // inserts particle in the cell neighbour list


    ~simpleCellNeighbourList();

    void sanityCheck();

    std::string reportCellNeighbourTopology() const;

    std::string reportParticleDistribution(const Real * positions, int iStart, int iEnd , int N) const;




    template<class V_t>
    Real twoBodyPotential(const double * positions, int iStart, int iEnd, int N,  const V_t & V, Real cutOff)
    {
        Real sum=0;
        if (dimensions == 3)
        {
                sum=evaluate_two_body_potential_triangular_unordered<3,V_t>
                (
                    positions, N , iStart , iEnd,
                    cellIndexPerParticle,
                    cells.data(), nParticlesPerCell.data()  ,
                    indexCellNeighbours.data(),
                    nNeighboursPerCell.data(),
                    V ,
                    cutOff
                );
        
        }

        return sum;
    }



    private:
 
    inline int index(int iX, int iY, int iZ) const {return iZ* nCells[0]*( nCells[1] + iY )  + iX  ;} // transform a 3d index in a 1d intex
    inline int index(int iX, int iY) const {return  nCells[0]* iY  + iX  ;} // transform a 2d index in a 1d intex
    inline int index(int iX) const {return  iX  ;} // transform a 1d index in a 1d intex

    inline int index(const std::array<int,3> & im) const {return index(im[0],im[1],im[2]);}


    int dimensions;
    realAllocator_t realAllocator;
    intAllocator_t intAllocator;

    // subdivision properties
    std::vector<Real> lBoxCell;
    std::vector<Real> lBox;
    std::vector<int> nCells;


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


};

}