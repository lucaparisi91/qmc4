
#ifndef PARTICLE_CONTAINER_H
#define PARTICLE_CONTAINER_H


#include "../traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "../toolsPimc.h"


namespace pimc
{

class cell 
{
    public:

    cell( size_t nMaxParticles=0); 

    void setCapacity(size_t nMaxParticles);


    public:

    
    const auto size() const  {return _nParticles;}

    const auto capacity() const  {return _nMaxParticles;}

    const auto & neighbours() const  {return _neighbours; }

    auto & neighbours()   {return _neighbours; }

    const auto & indices() const { return _particleIndex;};

    size_t push( int iParticle ) {

        if (_nParticles>=_particleIndex.size() )
        {
           setCapacity(_nParticles + nBuffer);
        }
        _particleIndex[_nParticles]=iParticle;

        _nParticles+=1;


        return _nParticles-1;

    }


    void swap( size_t ii , size_t jj  )
    {
        for(int d=0;d<getDimensions();d++)
        {
            std::swap(_particleIndex[ii],_particleIndex[jj]);
        }
    }

    void pop( ) { _nParticles-=1; }

    const auto & displacements() const {return _displacements;}
    
    auto & displacements()  {return _displacements;}

     void setCellIndex( const std::array<int,getDimensions()> & index) {
        _index[0]=index[0];
        _index[1]=index[1];
        _index[2]=index[2];
        

    }

    const auto & getCellIndex() const {return _index;}

    const auto & getParticleIndices() const {return _particleIndex;};

    
    private:

    size_t _nParticles;
    size_t _nMaxParticles;

    std::vector<int > _particleIndex;
    std::vector< cell * > _neighbours;
    Eigen::Tensor<Real,2 > _displacements;

    std::array<int,getDimensions()> _index;

    int nBuffer;


    /* std::Array<size_t,getDimensions()> index; */




};


inline auto pbc(Real x , Real L, Real LInverse)
{
    return ( x - std::round(x*LInverse )*L);

}


class simpleCellNeighbourList
{
    public:
    /*
    Naive data structure acceleration for short ranged ranged interactions.
    Subdivide the box in nCells[d] cell per dimension d.
    */

    simpleCellNeighbourList(std::array<size_t,getDimensions()> nCells, std::array<Real,getDimensions()> lBox   );

    auto capacity() const { return cellIndexPerParticle.size(); }

    auto size() const { return nParticles ; }

    void add( int  iParticle,   double x, double y, double z )
    {

        x=pbc(x,_lBox[0], _lBoxInverse[0]);
        y=pbc(y,_lBox[1], _lBoxInverse[1]);
        z=pbc(z,_lBox[2], _lBoxInverse[2]);


        auto i = std::floor( (x - _left[0] )/_delta[0] );
        auto j = std::floor( (y - _left[1]) /_delta[1] );
        auto k = std::floor( (z - _left[2]) /_delta[2] );

        auto I=index(i,j,k);

        cellIndexPerParticle[iParticle]=I;
        subIndexPerParticle[iParticle]=cells[I]->push(iParticle);
        isParticleRemoved[iParticle]=false;


        nParticles+=1;

    }


    void remove(  int iParticle)
    {
        if (not isParticleRemoved[iParticle] )
        {
            auto I=cellIndexPerParticle[iParticle];
            auto II=subIndexPerParticle[iParticle];
            cells[I]->swap( II, cells[I]->size()-1 );
        

        subIndexPerParticle[ cells[I]->indices()[II] ]=II;
        cells[I]->pop();
        isParticleRemoved[iParticle]=true;
        nParticles-=1;
        }

    }

    const auto & getParticleCell( int iParticle) const { auto I=cellIndexPerParticle[iParticle] ; return *cells[I];}

    bool isRemoved( int iParticle) const
    {
        return isParticleRemoved[iParticle];

    }
    
    void update(  size_t iParticle, Real x, Real y, Real z)
    {
        remove(iParticle);
        add(iParticle,x,y,z);
    };


    auto & getCell( size_t i, size_t j , size_t k ) { return *cells[index(i,j,k)];}

    const auto & getCell( size_t i, size_t j , size_t k ) const  { return *cells[index(i,j,k)];}


    ~simpleCellNeighbourList();

    void setCapacity( size_t N);

    const auto & getCells() { return cells;}

    const auto & getNCells() const {return _nCells;}

    private:
 
    inline int index(int iX, int iY, int iZ) const {return  _nCells[0]*( iZ* _nCells[1] + iY )  + iX  ;} // transform a 3d index in a 1d intex
    inline int index(int iX, int iY) const {return  _nCells[0]* iY  + iX  ;} // transform a 2d index in a 1d intex
    inline int index(int iX) const {return  iX  ;} // transform a 1d index in a 1d intex




    inline int index(const std::array<int,3> & im) const {return index(im[0],im[1],im[2]);}



    // subdivision properties

    std::array<Real,getDimensions() > _lBox;
    std::array<Real,getDimensions() > _left;
    std::array<Real,getDimensions() > _right;
    std::array<Real,getDimensions() > _delta;
    std::array<Real,getDimensions() > _lBoxInverse;
        
    std::array<size_t,getDimensions() > _nCells;

    // cell data
    std::vector< cell* > cells; 


    // per particle data
    std::vector< size_t > cellIndexPerParticle;
    std::vector< size_t > subIndexPerParticle;
    std::vector< bool > isParticleRemoved;
    
    int nParticles;

};

class maskTensor;

using mask_t=maskTensor;

class linkedCellParticles
{
    public:
    using linkedCell_t = pimc::simpleCellNeighbourList;

    linkedCellParticles( std::array<size_t,getDimensions()> nCells, std::array<Real,getDimensions() > lBox  );
    void setCapacity(size_t N, size_t M);
    void add( const Eigen::Tensor<Real,3> & data, const range_t &  timeRange,const  range_t & particleRange);
    void add( const Eigen::Tensor<Real,3> & data,const mask_t & mask, const range_t &  timeRange,const  range_t & particleRange);
    

    void update( const Eigen::Tensor<Real,3> & data, const range_t &  timeRange,const  range_t & particleRange);
     void update( const Eigen::Tensor<Real,3> & data, const mask_t & mask, const range_t &  timeRange,const  range_t & particleRange);

    void remove(  const range_t &  timeRange,const  range_t & particleRange);


    auto & operator[](size_t i ) {return *particles[i]; }
    const auto & operator[](size_t i ) const {return *particles[i]; }

    const auto & getNCells() const {return _nCells;}

    const auto & getBox() const {return _lBox;}

    const auto cellLength( int d) const {return _lBox[d]/_nCells[d];}

    private:


    std::vector< std::shared_ptr<linkedCell_t> > particles;
    std::array<size_t,getDimensions()> _nCells;
    std::array< Real ,getDimensions()> _lBox;

};


}


#endif