#include <vector>


namespace particleKernels
{

struct cellData
{
    virtual void updatePosition(int iSubCell,Real x,Real y , Real z){};

};



class indexCellData : public cellData
{   
    public:

    using index_t = int;


    indexCellData():
    _nParticles(0)
    {}

    void resize(index_t newSize)  {_particleIndices.resize(newSize) ;
    
    }

    size_t size() const noexcept {return _nParticles;}


    index_t add(index_t iParticle){
        if (_nParticles >=_particleIndices.size())
        {
            _particleIndices.resize(_nParticles + buffer);
        }
        _particleIndices[_nParticles]=iParticle;
        _nParticles++;
        return _nParticles-1;
    }

    auto nParticles() const  {return _nParticles;}

    void remove(index_t subCellIndex)
    {
        std::swap(_particleIndices[_nParticles-1] , _particleIndices[subCellIndex]);
        _nParticles--;
    }

    auto moveTo(int subCellIndex, indexCellData & cellDestination)
    {
        int iNewSubCell=cellDestination.add( getParticleIndex(subCellIndex ));
        remove(subCellIndex);
        return iNewSubCell;
    }  

    index_t  getParticleIndex(index_t i) const {
        return _particleIndices[i]; 
        } // get the i_th particle index

    private:  


    const int buffer = 10;
    std::vector<index_t> _particleIndices;
    index_t _nParticles; // number of particles containted in the box
};

class particleCellData3d : public cellData
{   
    public:

    using index_t = int;
    static constexpr int dims = 3;

    particleCellData3d():
    _nParticles(0)
    {}


    void resize(index_t newSize)  {
        _particleIndices.resize(newSize) ;
        for(int d=0;d<dims;d++)
        {
            _positions[d].resize(newSize);
        }
    }


    size_t size() const noexcept {return _nParticles;}

    index_t add(index_t iParticle, Real x,Real y, Real z){

        if (_nParticles >=_particleIndices.size())
        {
            resize(_nParticles + buffer);
        }
        _particleIndices[_nParticles]=iParticle;
        _positions[0][_nParticles]=x;
        _positions[1][_nParticles]=y;
        _positions[2][_nParticles]=z;

        _nParticles++;
        return _nParticles-1;
    }

    void updatePosition(int iSubCell,Real x,Real y , Real z) override
    {
        _positions[0][iSubCell]=x;
        _positions[1][iSubCell]=y;
        _positions[2][iSubCell]=z;
        
    };

    


    auto nParticles() const noexcept {return _nParticles;} 

    void remove(index_t subCellIndex)
    {
        std::swap(_particleIndices[_nParticles-1] , _particleIndices[subCellIndex]);


        for(int d=0;d<dims;d++)
        {
            std::swap(_positions[d][_nParticles-1] , _positions[d][subCellIndex]);
        }

        _nParticles--;
    }
  
    index_t inline   getParticleIndex(index_t i) const noexcept {
        return _particleIndices[i]; 
        } // get the i_th particle index
    

    Real inline   getParticlePosition(index_t i,int d) const noexcept {
        return _positions[d][i]; 
        } // get the i_th particle index

    
    auto moveTo(index_t subCellIndex, particleCellData3d & cellDestination)
    {
        int iNewSubCell=cellDestination.add( 
            getParticleIndex(subCellIndex ),
            getParticlePosition(subCellIndex,0 ),
            getParticlePosition(subCellIndex,1 ),
            getParticlePosition(subCellIndex,2 )
            );
        
        remove(subCellIndex);
        return iNewSubCell;
    }    
    private:

    const int buffer = 10;
    std::vector<index_t> _particleIndices;
    std::array<std::vector<Real>, dims> _positions;
    index_t _nParticles; // number of particles containted in the box
};

};

