
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

    const auto & positions() const {return _positions;};

    const auto size() const  {return _nParticles;}

    const auto capacity() const  {return _nMaxParticles;}

    const auto & neighbours() const  {return _neighbours; }

    auto & neighbours()   {return _neighbours; }

    const auto & indices() const { return _particleIndex;};


    size_t push( int iParticle, Real x , Real y , Real z ) {
        _positions(_nParticles,0)=x;
        _positions(_nParticles,1)=y;
        _positions(_nParticles,2)=z;
        _particleIndex[_nParticles]=iParticle;

        _nParticles+=1;

        return _nParticles-1;

    }


    void swap( size_t ii , size_t jj  )
    {
        for(int d=0;d<getDimensions();d++)
        {
            std::swap(_positions(ii,d),_positions(jj,d));
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


    
    private:

    size_t _nParticles;
    size_t _nMaxParticles;

    Eigen::Tensor<Real,2 > _positions;
    std::vector<int > _particleIndex;
    std::vector< cell * > _neighbours;
    Eigen::Tensor<Real,2 > _displacements;

    std::array<int,getDimensions()> _index;


    /* std::Array<size_t,getDimensions()> index; */




};






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
        auto i = std::floor( (x - _left[0] )/_delta[0] );
        auto j = std::floor( (y - _left[1]) /_delta[1] );
        auto k = std::floor( (z - _left[2]) /_delta[2] );

        auto I=index(i,j,k);

        cellIndexPerParticle[iParticle]=I;
        subIndexPerParticle[iParticle]=cells[I]->push(iParticle,x,y,z);

        nParticles+=1;

    }



    void remove(  int iParticle)
    {
        auto I=cellIndexPerParticle[iParticle];
        auto II=subIndexPerParticle[iParticle];
        cells[I]->swap( II, cells[I]->size()-1 );

        subIndexPerParticle[ cells[I]->indices()[II] ]=II;
        cells[I]->pop();

        nParticles-=1;

    }

    const auto & getParticleCell( int iParticle) { auto I=cellIndexPerParticle[iParticle] ; return *cells[I];}




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
    
    
        
    std::array<size_t,getDimensions() > _nCells;

    // cell data
    std::vector< cell* > cells; 


    // per particle data
    std::vector< size_t > cellIndexPerParticle;
    std::vector< size_t > subIndexPerParticle;
    

    
    int nParticles;
   
};

};
