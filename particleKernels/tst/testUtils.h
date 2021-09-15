#ifndef TEST_UTILS_H
#define TEST_UTILS_H

#include <random>

using Real = double;


class gaussianPotential
{
    public:
    gaussianPotential(Real sigma) : _alpha(1/(2*sigma*sigma)){}

    double operator()(double r) const {return exp(-_alpha*r*r);}
    double radialDerivative(double r) const {return -2*_alpha*r*exp(-_alpha*r*r);}

    auto alpha() {return _alpha;}
    private:
    double _alpha;
};



void generateRandomParticlesInBox(double* positions,int iStart,int iEnd,int N , int D, std::default_random_engine & generator, const double * lBox);
void generateRandomParticlesInBox(double* positions,int iStart,int iEnd,int t0 , int t1, int N , int D,int T, std::default_random_engine & generator, const double * lBox);

void setConstant(double* positions,int iStart,int iEnd,int t0 , int t1, int N , int D,int T, Real C);

void updateParticlesRandom(double* positions,int iStart,int iEnd,int N , int D, std::default_random_engine & generator, double delta=1);
void updateParticlesRandom(double* positions,int iStart,int iEnd,int t0,int t1,int N , int D,int T, std::default_random_engine & generator, double delta=1);

#endif