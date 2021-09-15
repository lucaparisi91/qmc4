
#include "particleData.h"
#ifndef TWO_BODY_DISTANCES_HPP
#define TWO_BODY_DISTANCES_HPP
#include "timers.h"

namespace particleKernels
{

struct distancePolicy
{


    inline auto static difference(const indexCellData & iCell, const indexCellData & jCell, int iSubCellIndex,int jSubCellIndex, const Real * particles, int d,int N,Real wrap) noexcept
    {
        auto iParticle = iCell.getParticleIndex(iSubCellIndex);
        auto jParticle = jCell.getParticleIndex(jSubCellIndex);
        return ( particles[iParticle + d*N ] - (particles[jParticle + d*N] 
                        + wrap 
                        )
        ) ;
    }

    inline auto static difference(const particleCellData3d & iCell, const particleCellData3d & jCell, int iSubCellIndex,int jSubCellIndex, const Real * particles, int d,int N,Real wrap) noexcept
    {

        auto iParticle = iCell.getParticleIndex(iSubCellIndex);
        auto jParticle = jCell.getParticleIndex(jSubCellIndex);
        //assert( iCell.getParticlePosition(iSubCellIndex,d) == particles[iParticle + d*N ] );
        //assert( jCell.getParticlePosition(jSubCellIndex,d) == particles[jParticle + d*N ] );


        //return ( particles[iParticle + d*N ] - (particles[jParticle + d*N] 
        //                + wrap 
        //                ));

        return ( iCell.getParticlePosition(iSubCellIndex,d) - (
             jCell.getParticlePosition(jSubCellIndex,d)
                        + wrap 
                        )
        ) ;
        

    } 


};

struct triangularLoopPolicy : public distancePolicy
{
    static bool accept(int i, int j, int iStart, int iEnd) noexcept
    {
        return not ( (j <= iEnd) and (j >= iStart) and (j>=i) );
    }
};

struct noParticleFilterPolicy : public distancePolicy
{
    static bool accept(int i, int j, int iStart, int iEnd) noexcept
    {
        return true;
    }
};


template<class particleData_t,class V_t,class policity_t = noParticleFilterPolicy>
Real twoBodyDistancesIsotropicReduction( const particleData_t & container , const V_t & op, const Real * particles , int iStart, int iEnd , int N   )
{
    Real sum2b=0;
   
    // loop over particles to be updated



{
    
    const auto & currentLattice = container.getLattice();
    const auto dimensions = currentLattice.dimensions();

#pragma omp parallel for reduction(+:sum2b) schedule(runtime)
    for (int iParticle=iStart;iParticle<=iEnd;iParticle++)
    {
        int iICell=container.cellIndex(iParticle);
        const auto & iCell = container[iICell];

        int iSubCell=container.subCellIndex(iParticle);

        //assert(iCell<currentLattice.extendedSize());
        //assert(iSubCell<container[iCell].size());

        for (int j=0;j<currentLattice. nCellsNeighbourhood() ;j++)
        {
            int iJCell=currentLattice.getNeighbour(iICell,j);
            //assert(jCell<currentLattice.extendedSize());
            const auto & jCell = container[iJCell];

            int nJcell=jCell.size();
            std::array<Real,3> warp {currentLattice.wrap(iJCell,0),currentLattice.wrap(iJCell,1),currentLattice.wrap(iJCell,2)};
            

#pragma omp simd reduction(+:sum2b)
            for (int jj=0;jj<nJcell;jj++)
            {
                int jParticle=jCell.getParticleIndex(jj);

                //assert(jParticle<N);

                if (policity_t::accept(iParticle,jParticle,iStart,iEnd) )
                {
                    Real r2=0;


                    for (int d=0;d<3;d++)
                    {
                        
                        Real diffd = policity_t::difference(iCell , jCell , iSubCell,jj, particles , d , N , warp[d] );
                         
                         
                        

                        //Real diffd=particles[iParticle+d*N] - particles[jj + d*N];

                        //diffd=utils::differencePBC(diffd,currentLattice.lengthBox()[d],1./currentLattice.lengthBox()[d] );
                        

                        
                        r2+=diffd*diffd;
                    }
                    sum2b+=op(std::sqrt(r2));
                }

            }

        }

    
    }
}
    
    return sum2b;
}

}
#endif