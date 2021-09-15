
#include "abstractFactory.h"
#include "estimators.h"
#include "estimatorCollection.h"
#include "correlationEstimator.h"

/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/

template<class T>
struct observableTraits;

typedef estimatorBase* (*estimatorCreatorFunc) ( const json_t & j);
typedef storer* (*storerCreatorFunc) ( const json_t & j);

realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j);
realScalarEstimator* createEstimatorFromOb(realScalarObservable * ob,const json_t & j);
realVectorEstimator* createEstimatorFromOb(realVectorObservable * ob,const json_t & j);

template<class observable_t>
estimatorBase * createEstimator(const json_t & j )
{
    
  
  observable_t * observable = new observable_t(j);
  
  auto est= createEstimatorFromOb(observable,j);
  
  return est;
}

template<class observable_t>
storer * createStorer(const json_t & j )
{
  observable_t * observable = new observable_t(j);
  
  auto storerPtr= new  typename estimatorTraits<typename observable_t::accumulator_t>::storer_t (observable,j);
  
  return storerPtr;
}

class estimatorFactory : public abstractFactory<estimatorBase,std::string, estimatorCreatorFunc>
{
public:
  using abstractFactory_t= abstractFactory<estimatorBase,std::string, estimatorCreatorFunc>;
  
  
  template<class ob_t >
  void registerObservable()
  {
    registerType( ob_t::name()  , & (createEstimator<ob_t> ) );
    
    knownObservableTypes[ob_t::name() ]=ob_t::kind();
    
  }
  
  std::vector<estimatorBase*> create(const json_t & j);
  
private:
  std::map<std::string,std::string> knownObservableTypes;
  
  
};



centerOfMassStorer* createSuperfluidFractionStorer(const json_t & j);


class storerFactory : public abstractFactory<storer,std::string, storerCreatorFunc>
{
public:
  using abstractFactory_t = abstractFactory<storer,std::string, storerCreatorFunc > ;

  
  template<class ob_t >
  void registerObservable()
  {
    registerType( ob_t::name()   , & (createStorer<ob_t> ) );
  }

 
 

  
  auto create(const json_t & j)
  {
    std::vector<storer*> storers;
    
    for (const auto & storeJson : j )
      {

	std::string id = storeJson["kind"].get<std::string>();
	
	
	if ( storeJson.find("recordSteps") !=storeJson.end() )
	  {
	    auto recordSteps = storeJson["recordSteps"].get<int >();
	    if ( recordSteps  > 0 )
	      {
		if ( (id != "forceEnergy") and (id !="energy") )
		  storers.push_back( abstractFactory_t::create(id,storeJson) );
		
	      }
	  }
	
      }

    centerOfMassStorer * cmptr = createSuperfluidFractionStorer(j);
    
    if ( cmptr != NULL )
      {
	
	storers.push_back(cmptr);
      }
    
    return storers;
  }
  
 
  
};




