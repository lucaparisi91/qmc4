
#include "particleDistancesList.h"
#include "particleContainer.h"
#include "geometryPMC.h"


namespace pimc
{

    twoBodyPairsList::twoBodyPairsList(size_t N,size_t M):
    filteredParticles( N,M+1)
    {
        nParticles.resize(M+1,0);
        currentParticleList.resize(N,0);
        nextParticleList.resize(N,0);
    };

    void twoBodyPairsList::buildParticleList( const linkedCellParticles & particles, size_t iParticle, range_t timeRange)
    {
        for (int t=timeRange[0];t<=timeRange[1]+1;t++)
        {
            int k=0;
            const auto & currentParticles=particles[t];

            const auto & cell=currentParticles.getParticleCell(iParticle);
            
            // loop on particles of the cell
            for( int jj=0 ; jj<cell.size();jj++ )
            {
                auto j = cell.indices()[jj];

                if (j!= iParticle)
                {
                    filteredParticles(k,t)=j;
                    k++;
                }
                

            }

            for ( auto neighbour : cell.neighbours() )
            {
                const auto & indices = neighbour->indices();
                for( int jj=0 ; jj<neighbour->size();jj++ )
                {
                    auto j = neighbour->indices()[jj];
                    filteredParticles(k,t)=j;
                    k++;

                }

            }


            nParticles[t]=k;
            
        }

        // sort particle indices in the whole time range


         for (int t=timeRange[0];t<=timeRange[1]+1;t++)
         {
             std::sort(  &filteredParticles(0,t) , &filteredParticles(nParticles[t],t) );
         }

         for (int t=timeRange[0];t<=timeRange[1];t++)
         {
             auto n = intersect( & filteredParticles(0,t), & filteredParticles(nParticles[t],t) , & filteredParticles(0,t+1), & filteredParticles(nParticles[t+1],t+1)  );
             nParticles[t]=n;
         }


    }



twoBodyDistancesBuffer::twoBodyDistancesBuffer(size_t N,size_t M) : pairList(N,M+1),
buffer(N,getDimensions(),M+1),
bufferOffset(N,getDimensions(),M+1)
{
    nDistances.resize(M+1,0);

}


void twoBodyDistancesBuffer::buildDistanceList( const  Eigen::Tensor<Real,3> & data, const linkedCellParticles &  cellList,  const geometry_t & geo , const range_t & timeRange, int iParticle)
{

    pairList.buildParticleList(cellList,iParticle,timeRange);

    const auto & nPairs = pairList.sizes();
    const auto & pairs = pairList.getPairList();

    int k=0;

    for(int t=timeRange[0];t<=timeRange[1];t++)
    {
        for( int ii=0;ii<nPairs[t];ii++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                buffer(k,d,t)=geo.difference( data(iParticle,d,t) - data(pairs(ii,t),d,t) ,d );
                bufferOffset(k,d,t)=data(iParticle,d,t+1) - data(pairs(ii,t),d,t+1);
            }
            
        }

    }




}







}