
using Real = double;


template<int dims,class V_t>
void addTwoBodyIsotropicForcesRectangular(
    const V_t & V,
    const Real * particles,
    Real * forces,
     int i1,int i2, // i-particle range
     int j1, int j2, // j-particle range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int D, // length of the second dimension of the array( usually equal to dims)
     int T, // total number of time slices,
     int NF, int DF , int TF, // timensions of the force array
     const Real * lBox,
     Real C // multiplication factor
     )
{
    /*
    Assumes no intersection between i-particle range and j-particle range 
    */
    Real sum2b=0;

    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for schedule(static) collapse(3)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<=j2;jParticle++)
            {
                Real r2=0;
                std::array<Real,3> diff;

                for(int d=0;d<dims;d++)
                {
                    diff[d]=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diff[d]= utils::differencePBC(diff[d],lBox[d],lBoxInverse[d]);
                    r2+= diff[d] * diff[d];
                }

                auto r = std::sqrt(r2);
                auto rInverse = 1./r;
                auto dVdr = V.radialDerivative(r);                
                
                for(int d=0;d<dims;d++)
                {
                    forces[iParticle + d*NF + t*N*DF ]+=dVdr*diff[d]*rInverse*C;
                    forces[jParticle + d*NF + t*NF*DF ]-=dVdr*diff[d]*rInverse*C;
                }
             }
        }
}

template<int dims,class V_t>
void addTwoBodyIsotropicForcesTriangular(
    const V_t & V,
    const Real * particles,
    Real * forces,
     int i1,int i2, // i-particle range
     int j1, // beginning of the i-range
     int t0 , int t1, // time range
     int N, // total number of particles
     int D, // max number of dimensions 
     int T, // total number of time slices,
     int NF, int DF, int TF , // dimensions of the 3d force array
     const Real * lBox,
     Real C // multiplication factor
     )
{
    /*
    Assumes no intersection between i-particle range and j-particle range 
    */
    Real sum2b=0;

    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }


    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                std::array<Real,3> diff;

                for(int d=0;d<dims;d++)
                {
                    diff[d]=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diff[d]= utils::differencePBC(diff[d],lBox[d],lBoxInverse[d]);
                    r2+= diff[d] * diff[d];
                }

                auto r = std::sqrt(r2);
                auto rInverse = 1./r;
                auto dVdr = V.radialDerivative(r);                
                
                for(int d=0;d<dims;d++)
                {
                    forces[iParticle + d*NF + t*NF*DF ]+=dVdr*diff[d]*rInverse*C;
                    forces[jParticle + d*NF + t*NF*DF ]-=dVdr*diff[d]*rInverse*C;
                }
             }
        }
}



template<int dims,class V_t>
Real evaluateTwoBodyRectangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1, int j2, // j-particle range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int D, // max number of dimensions
     int T, // total number of time slices
     const Real * lBox
     )
{
    Real sum2b=0;


    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<=j2;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                sum2b+=V(std::sqrt(r2));
        }
    }

return sum2b;
}


template<int dims,class V_t>
Real evaluateTwoBodyRectangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1, int j2, // j-particle range
     int t0 , int t1, // time range
     int N, // total number of particles 
     int D, // max number of dimensions
     int T, // total number of time slices
     const Real * lBox,
     const Real * mask
     )
{
    Real sum2b=0;

    std::array<Real,dims> lBoxInverse;
    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(3)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<=j2;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                Real ij_mask=  
                        0.5*( mask[ iParticle + (t+1)*N    ] * mask[ jParticle + (t+1)*N    ] + 
                        mask[ iParticle + (t)*N    ] * mask[ jParticle + (t)*N    ]

                         ) ;
                
                ij_mask = ((t == t0) or (t==t1)) ? std::min(ij_mask,0.5) : ij_mask ;



                sum2b+=V(std::sqrt(r2))*ij_mask;                

        }
    }

return sum2b;
}



template<int dims,class V_t>
Real evaluateTwoBodyTriangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1 ,// beginning of the i-particle set
     int t0 , int t1, // time range
     int N, // total number of particles
     int D , // maximum number of dimensions 
     int T, // total number of time slices
     const Real * lBox
     )
{
    Real sum2b=0;


    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(2)
    for (int t=t0;t<=t1;t++)
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                sum2b+=V(std::sqrt(r2));
        }
    }

return sum2b;
}


template<int dims,class V_t>
Real evaluateTwoBodyTriangular(
    const V_t & V,
    const Real * particles,
     int i1,int i2, // i-particle range
     int j1 ,// beginning of the i-particle set
     int t0 , int t1, // time range
     int N, // total number of particles
     int D , // maximum number of dimensions 
     int T, // total number of time slices
     const Real * lBox,
     const Real * mask
     )
{
    Real sum2b=0;


    std::array<Real,dims> lBoxInverse;

    for(int d=0;d<dims;d++)
    {
        lBoxInverse[d]=1./lBox[d];
    }

//#pragma omp parallel for reduction(+:sum2b) schedule(static) collapse(2)
    for (int t=t0;t<=t1;t++)
    {
        for (int iParticle=i1;iParticle<=i2;iParticle++)
        {
            for (int jParticle=j1;jParticle<iParticle;jParticle++)
            {
                Real r2=0;
                for(int d=0;d<dims;d++)
                {
                    Real diffd=( particles[iParticle + d*N + t*N*D ] - particles[jParticle + d*N + t*N*D] );
                    diffd= utils::differencePBC(diffd,lBox[d],lBoxInverse[d]);
                    r2+= diffd * diffd;
                }

                 Real ij_mask= 
                        0.5*( mask[ iParticle + (t+1)*N    ] * mask[ jParticle + (t+1)*N    ] + 
                        mask[ iParticle + (t)*N    ] * mask[ jParticle + (t)*N    ]

                         ) ;
                 ij_mask = ((t == t0) or (t==t1)) ? std::min(ij_mask,0.5) : ij_mask ;

                /* Real ij_mask=std::min( 
                        0.5*( mask[ iParticle + (t+1)*N    ] + mask[ iParticle + N*t    ] ), 
                          0.5*( mask[ jParticle + (t+1)*N    ] +  mask[ jParticle + N*t    ] ) 
                          ) ;
                 */
                sum2b+=V(std::sqrt(r2))*ij_mask;
                    
                    
                
            }
        }
    }


return sum2b;
}

template<int dims>
template<class V_t>
Real twoBodyPotential<dims>::operator()(const V_t & V, // two body potential
const Real * positions, // raw data in a contigous array of shape (N,dims, T)
int i1 , int i2, // particle range updated 
 int t0 , int t1  // time range updated
 // shape information of input data
 ) const
{

    int N = _dimensions[0];
    int D = _dimensions[1];
    int T = _dimensions[2];

    int iStartA = rangeA[0];
    int iStartB = rangeB[0];

    int iEndA = rangeA[1];
    int iEndB = rangeB[1];

    if (t0<0)
    {
        return 0;
    }

    Real sum=0;

    if ( not isTriangular )
    {
        int iA1= std::min(   std::max(i1,rangeA[0]) , rangeA[1] + 1 );
        int iA2= std::max(  std::min(i2,rangeA[1] ), rangeA[0]);

        int iB1= std::min( std::max(i1,rangeB[0]) , rangeB[0] );
        int iB2= std::max( std::min(i2,rangeB[1]) , rangeB[0] -1 );

        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iA1,iA2,rangeB[0],iB1-1,t0 , t1,N,D,T,_lBox.data());
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iA1,iA2,iB2+1,rangeB[1],t0 , t1,N,D,T,_lBox.data());
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,rangeA[0],iA1-1,t0 , t1,N,D,T,_lBox.data());
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,iA2+1,rangeA[1],t0 , t1,N,D,T,_lBox.data());


        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,iA1,iA2,t0 , t1,N,D,T,_lBox.data());
        
    }
    else
    {
        // interset the range of change particles [i1:i2] with rangeA
        i1=std::max(rangeA[0],i1);
        i2=std::min(rangeA[1],i2);

        sum+= evaluateTwoBodyTriangular<dims,V_t>(V,positions,i1,i2,iStartA,t0 , t1,N,D,T,_lBox.data());
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,i2+1,iEndA,i1,i2,t0 , t1,N,D,T,_lBox.data() );

    }

    return sum;

}

template<int dims>
template<class V_t>
Real twoBodyPotential<dims>::operator()(const V_t & V, // two body potential
const Real * positions, // raw data in a contigous array of shape (N,dims, T)
int i1 , int i2, // particle range updated 
 int t0 , int t1,  // time range updated
const Real * mask
 ) const
{
    
    int N = _dimensions[0];
    int D = _dimensions[1];
    int T = _dimensions[2];

    int iStartA = rangeA[0];
    int iStartB = rangeB[0];

    int iEndA = rangeA[1];
    int iEndB = rangeB[1];

    Real sum=0;

    if (t0<0)
    {
        return 0;
    }

    if ( not isTriangular )
    {
        int iA1= std::min(   std::max(i1,rangeA[0]) , rangeA[1] + 1 );
        int iA2= std::max(  std::min(i2,rangeA[1] ), rangeA[0]);

        int iB1= std::min( std::max(i1,rangeB[0]) , rangeB[0] );
        int iB2= std::max( std::min(i2,rangeB[1]) , rangeB[0] -1 );
        
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iA1,iA2,rangeB[0],iB1-1,t0 , t1,N,D,T,_lBox.data(),mask);
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iA1,iA2,iB2+1,rangeB[1],t0 , t1,N,D,T,_lBox.data(),mask);
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,rangeA[0],iA1-1,t0 , t1,N,D,T,_lBox.data(),mask);
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,iA2+1,rangeA[1],t0 , t1,N,D,T,_lBox.data(),mask);

        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,iB1,iB2,iA1,iA2,t0 , t1,N,D,T,_lBox.data(),mask);

    }
    else
    {
        // interset the range of change particles [i1:i2] with rangeA
        i1=std::max(rangeA[0],i1);
        i2=std::min(rangeA[1],i2);

        sum+= evaluateTwoBodyTriangular<dims,V_t>(V,positions,i1,i2,iStartA,t0 , t1,N,D,T,_lBox.data(),mask);
        sum+=evaluateTwoBodyRectangular<dims,V_t>(V,positions,i2+1,iEndA,i1,i2,t0 , t1,N,D,T,_lBox.data() ,mask);

    }

    return sum;

}


template<int dims>
template<class V_t>
void twoBodyPotential<dims>::addForce(
const V_t & V, // two body potential
const Real * positions, // raw data in a contigous array of shape (N,dims, T),
 Real * forces, // raw data in a contigous array of shape (N,dims, T) where to add the force for each particle
int i1 , int i2, // particle range updated 
 int t0 , int t1 , // time range updated,
 Real C // multiplication factor
 ) const
{
    int N= _dimensions[0];
    int D= _dimensions[1];
    int T= _dimensions[2];
    
    int NF= _forceDimensions[0];
    int DF= _forceDimensions[1];
    int TF= _forceDimensions[2];


    int iStartA = rangeA[0];
    int iStartB = rangeB[0];

    int iEndA = rangeA[1];
    int iEndB = rangeB[1];

   
    
    if ( not isTriangular )
    {
        // intersect the range of updated particles with rangeA and range B
        int iA1= std::min(   std::max(i1,rangeA[0]) , rangeA[1] + 1 );
        int iA2= std::max(  std::min(i2,rangeA[1] ), rangeA[0]);

        int iB1= std::min( std::max(i1,rangeB[0]) , rangeB[0] );
        int iB2= std::max( std::min(i2,rangeB[1]) , rangeB[0] -1 );

        // particle pairs AB form 5 rectangulars pair ranges. Compute forces in each particle pair subdomain
        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,iA1,iA2,rangeB[0],iB1-1,t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);
        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,iA1,iA2,iB2+1,rangeB[1],t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);
        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,iB1,iB2,rangeA[0],iA1-1,t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);
        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,iB1,iB2,iA2+1,rangeA[1],t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);

        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,iB1,iB2,iA1,iA2,t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);
        

    }
    else
    {
        i1=std::max(rangeA[0],i1);
        i2=std::min(rangeA[1],i2);

        addTwoBodyIsotropicForcesTriangular<dims,V_t>(V,positions,forces,i1,i2,iStartA,t0 , t1,N,D,T,NF,DF,TF,_lBox.data(),C);
        addTwoBodyIsotropicForcesRectangular<dims,V_t>(V,positions,forces,i2+1,iEndA,i1,i2,t0 , t1,N,D,T,NF,DF,TF,_lBox.data() ,C);
        
    }
}