
#include "../traits.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "../toolsPimc.h"



namespace pimc
{
    
    class simpleCellNeighbourList;

    using linkedCellList_t = simpleCellNeighbourList;

    class twoBodyPairsList
    {
        public:
        twoBodyPairsList(size_t N, size_t M);

        void buildParticleList( std::vector<linkedCellList_t> & particles , size_t iParticle, range_t timeRange); // generats a list of all particles interacting with the particle iParticle

        const auto & getPairList() {return filteredParticles; }
        const auto & sizes() {return nParticles; }


        private:


        Eigen::Tensor<int,2 > filteredParticles;
        std::vector<size_t> nParticles;

        std::vector< size_t > currentParticleList;
        std::vector< size_t > nextParticleList;



    };
}

