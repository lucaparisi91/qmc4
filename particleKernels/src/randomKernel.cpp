#include "randomKernel.h"

void randomSTLKernel::uniform(double * r , int n,double a,double b)
{
    std::uniform_real_distribution<double> udis(a,b);

    for (int i=0;i<n;i++)
    {
        r[i]=udis(randEngine);
    }

}

void randomSTLKernel::gauss(double * r , int n,double mean,double sigma)
{
    std::normal_distribution<double> ndis(mean,sigma);

    for (int i=0;i<n;i++)
    {
        r[i]=ndis(randEngine);
    }
}

#ifdef MKL

randomMKLKernel::randomMKLKernel(int seed) 
{
    vslNewStream( &stream, VSL_BRNG_MT19937, seed);
    
}

randomMKLKernel::~randomMKLKernel() 
{
    vslDeleteStream(&stream);
}

void randomMKLKernel::gauss(double * r , int n,double mean,double sigma)
{
    vdRngGaussian( VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, n, r, mean, sigma );
}

void randomMKLKernel::uniform(double * r , int n,double a,double b)
{
    vdRngUniform( VSL_RNG_METHOD_UNIFORM_STD, stream, n, r, a, b );
}

#endif