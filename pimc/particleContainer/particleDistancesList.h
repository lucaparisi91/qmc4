#ifndef PARTICLE_DISTANCES_H
#define PARTICLE_DISTANCES_H


#include "../traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "../toolsPimc.h"

namespace pimc
{
    class simpleCellNeighbourList;
    class linkedCellParticles;

    using linkedCellList_t = simpleCellNeighbourList;

    
    class twoBodyPairsList
    {
        public:
        twoBodyPairsList(size_t N, size_t M);

        void buildParticleList( const linkedCellParticles & particles , size_t iParticle, range_t timeRange); // generats a list of all particles interacting with the particle iParticle

        void buildParticleList( const linkedCellParticles & particles , const range_t & particleRange, range_t timeRange); // generate a list of all pair distances between particles with index in particle range

        const auto & sizes() {return nParticles; }

        const auto & getPairList() const {return filteredParticles;}



        private:

        Eigen::Tensor<int,3 > filteredParticles;
        

        std::vector<size_t> nParticles;

        std::vector< size_t > currentParticleList;
        std::vector< size_t > nextParticleList;

    };


    class twoBodyDistancesBuffer
    {
        public:

        twoBodyDistancesBuffer( size_t N,size_t M);
        void buildDistanceList( const  Eigen::Tensor<Real,3> & data, const linkedCellParticles &  cellList,  const geometry_t & geo , const range_t & timeRange, int iParticle);
         void buildDistanceList( const  Eigen::Tensor<Real,3> & data, const linkedCellParticles &  cellList,  const geometry_t & geo , const range_t & timeRange, const range_t & particleRange);

        const auto & getNDistances(size_t t) const { return nDistances[t];};
        const auto & getNDistances() const { return nDistances;};
        

        const auto & getBufferDistances() const { return buffer;};
        const auto & getBufferDistancesOffset() const { return bufferOffset;};
        
        const auto & getParticleList() const {return pairList;};
        const auto & getPairList() const {return _filteredParticles;};
        auto & getPairList()  {return _filteredParticles;};

    private:

    Eigen::Tensor<Real,3> buffer;
    Eigen::Tensor<Real,3> bufferOffset;

    std::vector<Real> nDistances;
    twoBodyPairsList pairList;
    Eigen::Tensor<int,3 > _filteredParticles;


};

}

#endif