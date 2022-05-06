
#include "particleDistancesList.h"
#include "particleContainer.h"
#include "geometryPMC.h"

namespace pimc
{
    twoBodyPairsList::twoBodyPairsList(size_t N,size_t M):
    filteredParticles( N*(N-1)/2,M+1,2)
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

            if (  not currentParticles.isRemoved(iParticle) )
            {
                        
                const auto & cell=currentParticles.getParticleCell(iParticle);
                
                // loop on particles of the cell
                for( int jj=0 ; jj<cell.size();jj++ )
                {
                    auto j = cell.indices()[jj];

                    if (j!= iParticle)
                    {
                        filteredParticles(k,t,1)=j;
                        k++;
                    }
                
                }


                for ( auto neighbour : cell.neighbours() )
                {
                    const auto & indices = neighbour->indices();
                    for( int jj=0 ; jj<neighbour->size();jj++ )
                    {
                        auto j = neighbour->indices()[jj];
                        filteredParticles(k,t,1)=j;
                        k++;

                    }

                }
            }

            nParticles[t]=k;

        }



    }


void twoBodyPairsList::buildParticleList( const linkedCellParticles & particles, const range_t & particleRange, range_t timeRange)
    {
        
        for (int t=timeRange[0];t<=timeRange[1]+1;t++)
        {
            int k=0;

            for(int iParticle=particleRange[0];iParticle<=particleRange[1];iParticle++)
            {
                const auto & currentParticles=particles[t];

                if (  not currentParticles.isRemoved(iParticle) )
                {
                        
                    const auto & cell=currentParticles.getParticleCell(iParticle);
                
                // loop on particles of the cell
                    for( int jj=0 ; jj<cell.size();jj++ )
                {
                    auto j = cell.indices()[jj];

                    if ( j < iParticle )
                    {
                        filteredParticles(k,t,0)=iParticle;
                        filteredParticles(k,t,1)=j;
                        k++;
                    }
                }


                for ( auto neighbour : cell.neighbours() )
                {
                    const auto & indices = neighbour->indices();
                    for( int jj=0 ; jj<neighbour->size();jj++ )
                    {
                        auto j = neighbour->indices()[jj];

                        if ( j< iParticle )
                            {
                                filteredParticles(k,t,0)=iParticle;
                                filteredParticles(k,t,1)=j;
                                k++;

                            }
                    }

                    }
                }

            
            }
            nParticles[t]=k;


        }

    }



twoBodyDistancesBuffer::twoBodyDistancesBuffer(size_t N,size_t M) : pairList(N,M+1),
buffer(N*(N-1)/2,getDimensions(),M+1),
bufferOffset(N*(N-1)/2,getDimensions(),M+1),
_filteredParticles(N*(N-1)/2,M+1,2)
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
                delta[d]=geo.difference( data(iParticle,d,t) - data(pairs(ii,t,1),d,t) ,d );
                deltaNext[d]=-(  - delta[d] - ( data(iParticle,d,t+1) - data(iParticle,d,t) ) +  data( pairs(ii,t,1) ,d,t+1) - data(  pairs(ii,t,1),d,t)); 
                r2+=delta[d]*delta[d];
                r2Next+=deltaNext[d]*deltaNext[d];
            }

            //std::cout << t << " " << pairs(ii,t) << " " << iParticle << " " << r2 << " " << r2Next << std::endl;
            
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

    }




}


void twoBodyDistancesBuffer::buildDistanceList( const  Eigen::Tensor<Real,3> & data, const linkedCellParticles &  cellList,  const geometry_t & geo , const range_t & timeRange, const range_t & particleRange)
{

    Real cutOff = cellList.cellLength(0);
    Real cutOffSquared = cutOff*cutOff;
    
    pairList.buildParticleList(cellList,particleRange,timeRange);

    
    
    
    const auto & nPairs = pairList.sizes();
    const auto & pairs = pairList.getPairList();
    auto & filteredList = getPairList();
    
    for(int t=timeRange[0];t<=timeRange[1];t++)
    {
        int k=0;

        for( int ii=0;ii<nPairs[t];ii++)
        {
            std::array<Real,getDimensions() > delta;
            std::array<Real,getDimensions() > deltaNext;
            Real r2=0;
            Real r2Next=0;

            auto iParticle = pairs(ii,t,0);
            auto jParticle= pairs(ii,t,1);

            for(int d=0;d<getDimensions();d++)
            {
                

                delta[d]=geo.difference( data(iParticle,d,t) - data(pairs(ii,t,1),d,t) ,d );
                deltaNext[d]=-(  - delta[d] - ( data(iParticle,d,t+1) - data(iParticle,d,t) ) +  data( pairs(ii,t,1) ,d,t+1) - data(  pairs(ii,t,1),d,t)); 
                r2+=delta[d]*delta[d];
                r2Next+=deltaNext[d]*deltaNext[d];
            }

            //std::cout << t << " " << pairs(ii,t) << " " << iParticle << " " << r2 << " " << r2Next << std::endl;
            
             if (  not((r2>=(cutOffSquared) ) or (r2Next>=cutOffSquared)
            ) )
            {
                for(int d=0;d<getDimensions();d++)
                {
                    buffer(k,d,t)=delta[d];
                    bufferOffset(k,d,t)=deltaNext[d];
                    filteredList(k,t,0)=iParticle;
                    filteredList(k,t,1)=jParticle;


                }
                k++;
            }

        }

        nDistances[t]=k;

    }




}





}