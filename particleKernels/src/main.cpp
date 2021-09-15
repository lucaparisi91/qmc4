using Real = double;

#include<cstdlib>
#include<iostream>
#include <numeric>
#include <immintrin.h>
#include "randomKernel.h"
#include "omp.h"

#define Time( name, command, NT) \
{ \
            double _ttime=0; \
            for (int _i = 0 ;_i< NT ; _i++ ) \
            { \
            double t1=omp_get_wtime(); \
            command ;\
            double t2=omp_get_wtime(); \
            _ttime+=t2-t1; \
            } \
            if ( omp_get_thread_num() == 0 )\
            {\
                std::cout << "Time\t " << name << "\t" <<  _ttime*1000 <<  "\t" << omp_get_max_threads() << std::endl;\
            }\
} 

Real harmonic_sum(const Real * data, int N )
{
    Real sum=0;

    for(int i=0;i<N;i++)
    {
        sum+=data[i]*data[i] + data[i+N]*data[i+N] + data[i+2*N]*data[i+2*N];
    }

    return sum;
}

#ifdef __AVX2__

Real harmonic_sum_avxd(const Real * data,int N)
{
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();

    for (int i=0;i<N/4;i++)
    {
        
        __m256d tmp1 = _mm256_load_pd(&data[i*4]);
        sum1 = _mm256_fmadd_pd(tmp1,tmp1,sum1);

        __m256d tmp2 = _mm256_load_pd(&data[i*4+N]);
        sum2 = _mm256_fmadd_pd(tmp2,tmp2,sum2);

        __m256d tmp3 = _mm256_load_pd(&data[i*4+2*N]);
        sum3 = _mm256_fmadd_pd(tmp3,tmp3,sum3);

    }

    sum1=_mm256_add_pd(sum1,sum2);
    sum1=_mm256_add_pd(sum1,sum3);

    double * res = (double*) &sum1;

    return res[0] + res[1] + res[2] + res[3];


}

Real harmonic_sum_avxd2(const Real * data,int N)
{
    __m256d sum1 = _mm256_setzero_pd();
    __m256d sum2 = _mm256_setzero_pd();
    __m256d sum3 = _mm256_setzero_pd();

    for (int i=0;i<N;i=i+4 )
    {
        
        __m256d tmp1 = _mm256_load_pd(&data[i]);
        sum1 = _mm256_fmadd_pd(tmp1,tmp1,sum1);
    }
    
    for (int i=N;i<2*N;i=i+4 )
    {
        
        __m256d tmp1 = _mm256_load_pd(&data[i]);
        sum2 = _mm256_fmadd_pd(tmp1,tmp1,sum2);
    }

    for (int i=2*N;i<3*N;i=i+4 )
    {
        
        __m256d tmp1 = _mm256_load_pd(&data[i]);
        sum3 = _mm256_fmadd_pd(tmp1,tmp1,sum3);
    }


    sum1=_mm256_add_pd(sum1,sum2);
    sum1=_mm256_add_pd(sum1,sum3);

    double * res = (double*) &sum1;

    return res[0] + res[1] + res[2] + res[3];



}

#endif





Real pair_potential_harmonic_sum( const Real *  data,int N )
{
    Real sum=0;

    for(int i=0;i<N;i++)
    {
        for(int j=0;j<i;j++)
        {
            
            sum+= (data[i] - data[j])*(data[i]-data[j]) + 
                 ( data[i+N]-data[j+N] ) * (data[i+N]-data[j+N]) +
                 ( data[i+2*N]-data[j+2*N] ) * (data[i+2*N]-data[j+2*N]);

        }
    }

    return sum;
}


class harmonic_potential
{
    public :
    harmonic_potential() {}
    inline double operator()(double x, double y , double z) const {return x*x + y*y + z*z;}
};

template<class V_t>
void pair_potential_sum( const Real *  data,int N , const V_t & V ,  double & sum)
{


    for(int i=0;i<N;i++)
        for(int j=0;j<i;j++)
        {
            sum+= V(
                 (data[i] - data[j]) ,
                 ( data[i+N]-data[j+N] ) ,
                 ( data[i+2*N]-data[j+2*N] )  
                );
        }
}


template<class V_t>
void potential_one_body_sum( const Real *  data,int N , const V_t & V ,  double & sum)
{
    for(int i=0;i<N;i++)
        {

            sum+= V(
                 (data[i] ) ,
                 ( data[i+N] ) ,
                 ( data[i+2*N] )  
                );
        }
}




void computeDriftForceHarmonic_naive(const Real * data1, const Real * data2, Real * forces1, Real * forces2, int N1,int N2)
{
    for(int i=0;i<N1;i++)
    {
        for(int j=0;j<N2;j++)
        {
            Real tmp = 0.5*(data1[i] - data2[j] );
            forces1[i]+=tmp;
            forces2[j]-= tmp;
        }

    }

}

void computeDriftForceHarmonic_opt(const Real * __restrict__ data1, const Real * __restrict__  data2, Real * __restrict__  forces1, Real * __restrict__  forces2, int N1,int N2)
{
    for(int i=0;i<N1;i++)
    {
        for(int j=0;j<N2;j++)
        {
            Real tmp = 0.5*(data1[i] - data2[j] );
            forces1[i]+=tmp;
            forces2[j]-= tmp;
        }

    }

}

int main(int argc,char ** argv)
{
    const long int N=32*1000;
    const int D=3;

    Real * data = (Real *)aligned_alloc( 32, N*D*sizeof(Real));
    Real * forces = (Real *)aligned_alloc( 32, N*D*sizeof(Real));
    Real * forces2 = (Real *)aligned_alloc( 32, N*D*sizeof(Real));
    Real * data2 = (Real *)aligned_alloc( 32, N*D*sizeof(Real));

    srand(115);

    for(int i=0;i<D*N;i++)
    {
        data[i]=rand()*1./RAND_MAX;
        data2[i]=rand()*1./RAND_MAX;
    }

/*     Real sum=0;
    Time ("harmonic_sum_avx" ,  sum = harmonic_sum( data, N ) )
    std::cout << "Sum: " << sum << std::endl; */
    /*

    sum=0;
    Time ("harmonic_sum_naive" ,  sum = harmonic_sum( data, N ) )
    std::cout << "Sum: " << sum << std::endl;

    sum=0;
    Time ("harmonic_sum_naive_avx2" ,  sum = harmonic_sum_avxd2( data, N ) )
    std::cout << "Sum: " << sum << std::endl;

    */
/* 
    
   Time ("force_opt" ,  computeDriftForceHarmonic_opt(data,data2,forces,forces2,N,N) )
    //Time ("force_naive" ,  computeDriftForceHarmonic_naive(data,data2,forces,forces2,N,N) )

    Real sum=0;
    for(int i=0;i<D*N;i++)
    {
        sum+=forces[i] + forces2[i];
    }

    std::cout << sum << std::endl; */

  /*   randomSTLKernel randKernelSTL(15);
    randomMKLKernel randKernelMKL(15);

    Time( "rand_STL" ,randKernelSTL.gauss(data,N*D) );
    Time( "rand_MKL" ,randKernelMKL.gauss(data2,N*D) );


    Real sum=0,sum2=0;

    for(int i=0;i<D*N;i++)
    {
       sum+=data[i]*data[i];
       sum2+=data2[i]*data2[i];
    }

    std::cout << sum/(N*D) << std::endl;
    std::cout << sum2/(N*D) << std::endl; */

    Real sum,sum2b=0,sum1b=0;

    harmonic_potential V;


    //Time ( "pair_harmonic_sum", sum+=pair_potential_harmonic_sum( data, N );, 1 )
    

    Time ( "pair_harmonic_sum_opt", pair_potential_sum( data2, N , V ,sum2b); , 10)
    
    Time ( "potential_one_body_sum", potential_one_body_sum( data, N , V ,sum1b); , 10)


    //std::cout << "sum: " << sum/(N*(N-1)*0.5) <<  std::endl;
    std::cout << "sum2: " << sum2b/(N*(N-1)*0.5) <<  std::endl;
    std::cout << "sum1b: " << sum1b/N <<  std::endl;


    free( data);





}