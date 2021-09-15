

template<int dimensions,class V_t>
double evaluate_two_body_potential_triangular_unordered(const double * positions, int N, int iStart, int iEnd,
const int * cellIndexPerParticle,
 const int * const * particlesPerCellIndex  , int * nParticlesPerCellIndex,
  const int * const * indexNeighborCellsPerCell , 
  const int * nNeighboursPerCellIndex,
const V_t & V, Real cutOff)
{
    double sum=0;

    for (int i=iStart;i<=iEnd ; i++ )
    {
        int iCell=cellIndexPerParticle[i];


        for (int jjCell=0;jjCell<nNeighboursPerCellIndex[iCell] ; jjCell++  )
        {
            int jCell=indexNeighborCellsPerCell[iCell][jjCell];


            for (int iParticle=0; iParticle<nParticlesPerCellIndex[jCell];iParticle++)
            {
                int j =particlesPerCellIndex[jCell][iParticle];

                if (j<i)
                {
                    double distanceSquared=0;               

                    for (int d=0; d< dimensions ; d++)
                        distanceSquared+=  (positions[i + N*d] - positions[j + N*d])*(positions[i + N*d] - positions[j + N*d] );
                    
                    double distance=std::sqrt(distanceSquared);

                    if (distance <= cutOff)
                    {
                        sum+=V(distance);
                    }

                    

                    
                }

            
            }
        }
    }

    return sum;
}