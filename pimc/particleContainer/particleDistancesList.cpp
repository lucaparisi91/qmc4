
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

        //std::cout << "..not intersected neighbours..." << std::endl;

        // sort particle indices in the whole time range
        /*  for (int t=timeRange[0];t<=timeRange[1]+1;t++)
         {
             std::cout << t << " ";

             for(int ii=0;ii<nParticles[t];ii++)
             {
                 std::cout  << filteredParticles(ii,t) << " ";
             }

             std::cout << std::endl;

         } */

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

    Real cutOff = cellList.cellLength(0);
    Real cutOffSquared = cutOff*cutOff;

    pairList.buildParticleList(cellList,iParticle,timeRange);

    const auto & nPairs = pairList.sizes();
    const auto & pairs = pairList.getPairList();

    for(int t=timeRange[0];t<=timeRange[1];t++)
    {
        int k=0;
        for( int ii=0;ii<nPairs[t];ii++)
        {
            std::array<Real,getDimensions() > delta;
            std::array<Real,getDimensions() > deltaNext;
            Real r2=0;
            Real r2Next=0;
            
            for(int d=0;d<getDimensions();d++)
            {
                delta[d]=geo.difference( data(iParticle,d,t) - data(pairs(ii,t),d,t) ,d );
                deltaNext[d]=geo.difference(data(iParticle,d,t+1) - data(pairs(ii,t),d,t+1),d);  
                r2+=delta[d]*delta[d];
                r2Next+=deltaNext[d]*deltaNext[d];
            }

             if (  not((r2>=(cutOffSquared) ) or (r2Next>=cutOffSquared)
            ) )
            {
                for(int d=0;d<getDimensions();d++)
                {
                    buffer(k,d,t)=delta[d];
                    bufferOffset(k,d,t)=deltaNext[d];
                }
                k++;
            }

        }

        nDistances[t]=k;

        
        //std::cout<<std::endl;


    }




}







}