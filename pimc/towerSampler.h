
#ifndef TOWERSAMPLER_H
#define TOWERSAMPLER_H


#include "tools.h"

namespace pimc
{
class towerSampler
{
    public:
    
    towerSampler(int max_num_weights) : _accWeights(max_num_weights,0.),uniformRealNumber(0,1), _totWeight(0),iCurrentWeight(0) {}

    void resize(int max_num_weights)  {_accWeights.resize(max_num_weights,0.);} 
    


    towerSampler() : towerSampler(0) {} 
   

    virtual void reset() {_totWeight=0;iCurrentWeight=0;_accWeights.resize(0);}

    
    
    void accumulateWeight(Real weight) {
        if ( iCurrentWeight >= _accWeights.size() )
        {
            _accWeights.resize(iCurrentWeight+1);    
        }
        _totWeight+=weight;
        _accWeights[iCurrentWeight]=_totWeight;
        iCurrentWeight++;
    }

    int sample(randomGenerator_t & randG)
    {
        return sample(_accWeights,_totWeight,randG);
    }

    int sample(std::vector<Real> & accWeights, Real totWeight, randomGenerator_t & randG)
    {
        Real randomWeight=uniformRealNumber(randG) * totWeight;

         if (iCurrentWeight== 0)
        {
            throw invalidState("Nothing to sample.");
        }
        int i=0;
        for(i=0; accWeights[i]<randomWeight;i++) {}

       
        return i;
    }

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
    std::vector<Real> _accWeights;
    Real _totWeight;
    int iCurrentWeight;

};

}

#endif