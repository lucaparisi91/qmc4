
#include "geometryPMC.h"
#include "tools.h"

namespace pimc
{
    
    void geometryPBC_PIMC::updateEqualTimeDifferences( Eigen::Tensor<Real,3>  & differences, const  Eigen::Tensor<Real, 3> & tn, std::array<int,2> timeRange,  int j,  std::array<int,2> particleRangeB )
    {
        for (size_t t = timeRange[0] ; t<= timeRange[1] ;t++  )
        {
            for ( size_t i=particleRangeB[0]  ;i<=particleRangeB[1];i++ )
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        differences(  i  ,d , t  )  = difference(      tn(i,d, t ) - tn(j,d,t) , d);

                    }
                    
                }
        }

    }

void geometryPBC_PIMC::updateSpringDifferences( Eigen::Tensor<Real,3> & differences , const  Eigen::Tensor<Real, 3> & tn, std::array<int,2> timeRange, std::array<int,2> particleRange )
    {
        size_t N = particleRange[1] - particleRange[0]+ 1;
        size_t T = timeRange[1] - timeRange[0] + 1;

        for (size_t t = timeRange[0] ; t<= timeRange[1] ;  t++ )
        {
            for ( size_t i=particleRange[0]  ;i<=particleRange[1];i++)
                {
                    for (int d=0;d<getDimensions();d++)
                    {
                        differences( t, d , i )  = 
                         difference(      tn(i,d, (t+1) ) - tn(i,d,t) , d);

                    }
                    
                }
                
                
        }
        
    }

    
}






// old code





/* 
geometryPBC_PMC::springDistances_t geometryPBC_PMC::springDistances(  const std::vector<states_t> & particles ) const
{
    int timeRangess = particles.size();

    springDistances_t differences;
    int nSets=particles[0].size();

    differences.resize(nSets);


    for (int iSet =0 ; iSet< nSets;iSet++)
    {
        int N=particles[0][iSet].rows();
        differences[iSet].resize(N);
        for(int n=0;n<N;n++)
        {
            differences[iSet][n].resize(timeRangess , getDimensions() );
        }
        for (int i = 0 ;i<timeRangess;i++)
        {
            const auto & X = particles[i][iSet];
            const auto & XN = particles[ (i+1)%timeRangess][iSet];

            for(int n=0;n<N;n++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    differences[iSet][n](i,d) =  difference( X(n,d) - XN(n,d) ,d) ;
                }
                
            }
                
            
        }
    }

    return differences;
}

void geometryPBC_PMC::springDistances( springDistances_t & distances, const std::vector<states_t> & particles, int startSlice , int endSlice, int iSet, int iParticle ) const
{
     auto & chain = distances[iSet][iParticle];
    int timeRangess = particles.size();

    for (int i=startSlice ; i< endSlice;i++)
    {
        for (int d=0;d<getDimensions();d++)
        {
        chain(i,d) = difference( particles[(i+1)%timeRangess][iSet](iParticle,d) - particles[i][iSet](iParticle,d)  ,d);
        }
    }

}


void geometryPBC_PMC::differencesTwoBody( difference_t &diffs, const geometryPBC::particles_t & particleData , int i, int j_start, int j_end) const 
	{
        // copute the diffeerences between particles i and particles j
	    const int N = particleData.rows();

        constexpr int D=difference_t::ColsAtCompileTime;

		// loop on a row of differences matrix
 		for(int j=j_start;j<=j_end;j++)
 		{
 			for (int d=0;d<D;d++)
			 {
 				diffs(j,d)=difference( particleData(i,d) - particleData(j,d) , d);
 				j++;
 			}	
		 }

	}


     */