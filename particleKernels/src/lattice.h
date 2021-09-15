#ifndef LATTICE_H
#define LATTICE_H


#include <vector>
#include <memory>
#include <cmath>
#include <cassert>

namespace particleKernels
{

using Real = double;

class lattice
{

    public:
    /*
    Lattice with with dimensions 1,2 or 3 and nGhosts[] ghost cells in each direction
    */

    using index_t = int;


    lattice(std::array<size_t,3>  shape, std::array<Real,3> lowerEdge, std::array<Real,3 > higherEdge,std::array<size_t,3> nGhosts={1,1,1}); // constructs a 3d lattice with box length lBox center on the origin

    inline int index(int i , int j , int k) const { return _extendedShape[0]*(  _extendedShape[1]*k + j   )  + i ; } // fuse a 3d index in a 1d index ( fortran ordering )

    inline int index(int i, int j) const {return _extendedShape[0]*j + i;}
    inline int index(int i) const {return i;}

    const auto & lengthBox() const {return _lBox;}
    const auto & inverseLengthBox() const {return _lBoxInverse;}
    const auto & lowerEdge() const {return _lowerEdge;}

    // map d indices to a single index
    inline index_t index(const std::array<int,3> & i) const {return index(i[0],i[1],i[2]);}

    inline index_t index(const std::array<int,2> & i) const {return index(i[0],i[1]);} 

    inline index_t index(const std::array<int,1> & i) const {return index(i[0]);}

    std::array<index_t,3> index3d(index_t iCell) const 
    {
        int k = iCell/(_extendedShape[0]*_extendedShape[1]);
        iCell=iCell - k*_extendedShape[0]*_extendedShape[1];
        int j = iCell/(_extendedShape[0]);
        int i = iCell - j*_extendedShape[0];
        return std::array<index_t,3>{i,j,k};
    }


    // map particle position to cell index
    inline index_t index(double x, int d) const noexcept{
         return std::floor( (  (x- _lowerEdge[d])*_lBoxInverse[d]  ) * _shape[d]) + _nGhosts[d];
         //assert( (i<=highIndex(d) ) & (i>=lowIndex(d)) );return i;
          } // x is the d coordinate of the dims-dimensional position vector (x,y,z)
    

    inline index_t index(double x, double y, double z) const { 
        return index( index(x,0), index(y,1),index(z,2)  );
    }
    inline index_t index(double x, double y) const {return index( index(x,0), index(y,1)  );}


    int dimensions() const noexcept{return dims;}
    const auto & shape() const noexcept{return _shape;}


    const auto & extendedShape() const noexcept {return _extendedShape;} // including ghost cells

    const auto & nGhosts() const {return _nGhosts;}



    size_t size() const noexcept {return _size;} // returns the total number of valid cells

    auto extendedSize() const noexcept {return _extendedSize;}
    
    int nCellsNeighbourhood() const noexcept {return _nNeighboursPerCell;}


    int getNeighbour(index_t iCell, index_t iNeighboorhood ) const noexcept  {return _neighbourOffsets[iNeighboorhood] + iCell; } // includes the cell i from some non zero offset


    inline int lowIndex(int d) const noexcept {return _nGhosts[d];}
    inline int highIndex(int d) const noexcept {return _shape[d] - 1 + _nGhosts[d];}



    inline Real wrap( int iCell,int d) const noexcept { return _wrap[d][iCell];} // return -1/0/1 * lBox[d] depending if the cell is on the low face, bulk , high face in direction d

    void checkNeighbourIndexing();
    

    private:

    void initialize();
    void buildIndexOffsets();

    void buildBC();
  
    const int dims;
    const bool periodic;
    const int _nNeighboursPerCell;
    
    // domain data

    std::array<Real,3> _lBoxInverse;
    std::array<Real,3 > _lBox;
    std::array<Real,3 > _lowerEdge;
    std::array<Real,3 > _higherEdge;
    size_t _size;
    size_t _extendedSize;



    // topological data
    std::array<size_t, 3> _shape;
    std::array<size_t, 3> _nGhosts;
    std::array<size_t,3> _extendedShape;

    std::array<index_t, 27> _neighbourOffsets;

    std::array<std::vector<Real> ,3    > _wrap;
};
}
#endif