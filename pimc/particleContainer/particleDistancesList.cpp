
#include "particleDistancesList.h"
namespace pimc
{

    twoBodyPairsList::twoBodyPairsList(size_t N,size_t M):
    filteredParticles( N,M)
    {
        nParticles.resize(M,0);
        currentParticleList.resize(N,0);
        nextParticleList.resize(N,0);
    };
    



}