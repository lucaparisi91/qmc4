#ifndef TOOLSPIMC_H
#define TOOLSPIMC_H

#include "traits.h"
#include "tools.h"
#include <random>
#include <array>

namespace pimc{



#if DIMENSIONS==1
    #define TRUNCATE_D(a,b,c) a
#endif

#if DIMENSIONS==2
    #define TRUNCATE_D(a,b,c) a,b
#endif

#if DIMENSIONS==3
    #define TRUNCATE_D(a,b,c) a,b,c
#endif



    using json_t = nlohmann::json;

    int getTotParticlesFromJson(const json_t & j);

    
    using Real = double;

    enum periodicity {periodic = 1 , open = 0};

    std::array<std::array<int,2>, 2> splitPeriodicTimeSlice(const std::array<int,2> & timeSlice, int nBeads);
    
    Real freeParticleLogProbability(std::array<Real,getDimensions()> & delta,Real tau,Real mass=1);
    


    Real average(const std::vector<Real> & observables);
    
}

#endif