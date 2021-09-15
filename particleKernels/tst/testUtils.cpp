#include "testUtils.h"

void generateRandomParticlesInBox(double* positions,int iStart,int iEnd,int t0 , int t1, int N , int D,int T, std::default_random_engine & generator, const double * lBox)
{
    // fill positions with uniform random particle positions [-lBox/2,lBox/2]

    std::uniform_real_distribution<double> unifRandom(-0.5,0.5);
    for(int t=t0;t<=t1;t++)
    for(int i=iStart;i<=iEnd;i++)
    {
        for(int d=0;d<D;d++)
        {
            positions[i+N*d + D*N*t]=unifRandom(generator)*lBox[d] ;
        }
    }
}



void generateRandomParticlesInBox(double* positions,int iStart,int iEnd,int N , int D, std::default_random_engine & generator, const double * lBox)
{
    // fill positions with uniform random particle positions [-lBox/2,lBox/2]

    std::uniform_real_distribution<double> unifRandom(-0.5,0.5);
    for(int i=iStart;i<=iEnd;i++)
    {
        for(int d=0;d<D;d++)
        {
            positions[i+N*d]=unifRandom(generator)*lBox[d] ;
        }
    }
}


void setConstant(double* positions,int iStart,int iEnd,int t0 , int t1, int N , int D,int T, Real C)
{
    //set with a constat C positions[iStart:iEnd,0:D-1,t0:t1] 
    std::uniform_real_distribution<double> unifRandom(-0.5,0.5);
    for(int t=t0;t<=t1;t++)
    for(int i=iStart;i<=iEnd;i++)
    {
        for(int d=0;d<D;d++)
        {
            positions[i+N*d + D*N*t]=C;
        }
    }
}



void updateParticlesRandom(double* positions,int iStart,int iEnd,int N , int D, std::default_random_engine & generator, double delta)
{
    std::uniform_real_distribution<double> unifRandom(0.0,1.0);
    
    for(int i=iStart;i<=iEnd;i++)
    {
        for(int d=0;d<D;d++)
        {
            positions[i+N*d]=+unifRandom(generator)*delta;
        }
    }
}


void updateParticlesRandom(double* positions,int iStart,int iEnd, int t0, int t1, int N , int D, int T,  std::default_random_engine & generator, double delta)
{
    std::uniform_real_distribution<double> unifRandom(0.0,1.0);
    
    for(int t=t0;t<=t1;t++)   
        for(int i=iStart;i<=iEnd;i++)
            {
                for(int d=0;d<D;d++)
                {
                    positions[i+N*d + t*N*D]=+unifRandom(generator)*delta;
                }
            }
}