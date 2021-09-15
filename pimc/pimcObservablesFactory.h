#ifndef PIMC_OBSERVABLES_FACTORY_H
#define PIMC_OBSERVABLES_FACTORY_H

#include "pimcObservables.h"

namespace pimc
{

template<class T>
std::shared_ptr<observable> __createObservable(const json_t & j)
 {
    auto estimator=std::make_shared<T>(j);

    auto ob = std::make_shared<typename T::observable_t>(estimator,j);
    
    return ob;

 }

class pimcObservablesFactory
{
    public:

    pimcObservablesFactory(int nBeadsMax,std::vector<int> nParticles);
    
    std::shared_ptr<observable> createObservable(const json_t & jOuter);


    std::vector<std::shared_ptr<observable> >  createObservables(const json_t & j);



    template<class T>
    void registerObservable(const std::string & key)
    {
        creatorMap[key]= &__createObservable<T> ; 
    }

    auto getEnergyObservable(){return eO;}
    private:

    int _nBeadsMax;
    std::vector<int> _nParticles;
    int _nChains;

    typedef std::shared_ptr<observable> (*observableCreatorFunc) ( const json_t & j);

using creatorMap_t = std::map<std::string,observableCreatorFunc>;
    creatorMap_t creatorMap;


    std::shared_ptr<scalarObservable> eO=nullptr;
};

}

#endif