#include "pimcObservablesFactory.h"

namespace pimc
{

 pimcObservablesFactory::pimcObservablesFactory(int nBeadsMax,std::vector<int> nParticles):
 _nBeadsMax(nBeadsMax),
 _nParticles(  nParticles )
 {
    _nChains = std::accumulate(_nParticles.begin(),_nParticles.end(),0);
 }

 
 std::shared_ptr<observable>  pimcObservablesFactory::createObservable(const json_t & jOuter)
    {
        // creates a copy of the json input and supplement with additional essential information
        json_t j(jOuter);
        j["nChains"]=_nChains;
        j["nBeads"]=_nBeadsMax;
        j["particles"]=_nParticles;


        std::string key=j["kind"].get<std::string>();

        creatorMap_t::const_iterator i;
        i=creatorMap.find(key);


        if (i!= creatorMap.end())
        {
	    	auto ob=(i->second)( j  );
            if (key=="thermalEnergy" )
            {
                eO=std::dynamic_pointer_cast<scalarObservable> (ob);
            }
            return ob;
        }
         else
        {
	    throw factoryIdNotRecorded(key);
        }

    }

std::vector<std::shared_ptr<observable> > pimcObservablesFactory::createObservables(const json_t & jOuter)
{
    std::vector<std::shared_ptr<observable> > observables;
    
    for (auto jO : jOuter)
    {

        observables.push_back(createObservable(jO));

    }

    return observables;
}

    

}