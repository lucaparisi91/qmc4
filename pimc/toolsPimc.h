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
    
    using range_t = std::array<int,2>;

    range_t intersectRanges(const range_t & rangeA,const range_t & rangeB );


    template<class IT1_t , class IT2_t>
    auto intersect( IT1_t left1, IT1_t right1,  IT2_t left2, IT2_t right2 )
    {
        auto head = left1;
        int n=0;

        while (left1 < right1 )
        {
            while ( *left2 < *left1 )
            {
                left2++;
                if (left2 >=right2 ) return n;
            }

            if (*left1 == *left2)
            {
                *head = *left1;
                head++;
                n++;
            }

            left1++;


        }
        
        return n;


    }
    
}

#endif