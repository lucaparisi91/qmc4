#include <random>
#include <random>



class randomSTLKernel
{   
    public:

    randomSTLKernel(int seed) : randEngine(seed)  {}
    void uniform(double * r , int n, double a=0 , double b=1); // generaters n uniform random numbers between 0 and 1 and save at memory location r
    void gauss(double * r, int n, double mean=0,double sigma=1); // generates n gaussian random number and saves at memory location r


    private:

    std::ranlux24_base randEngine;

};


#ifdef MKL
#include "mkl/mkl_vsl.h"


class randomMKLKernel
{
    public:

    randomMKLKernel(int seed) ;
    void uniform(double * r , int n,double a=0,double b=1.); // generaters n uniform random numbers between 0 and 1 and save at memory location r
    void gauss(double * r , int n,double mean=0.,double sigma=1.); // generaters n uniform random numbers between 0 and 1 and save at memory location r
    ~randomMKLKernel();


    private:


     VSLStreamStatePtr stream;



};

#endif
