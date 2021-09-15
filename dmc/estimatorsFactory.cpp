#include "estimatorsFactory.h"
#include "correlationEstimator.h"


realHistogramEstimator* createEstimatorFromOb(realHistogramObservable * ob,const json_t & j)
{
  return new realHistogramEstimator(ob,j);
}

realScalarEstimator* createEstimatorFromOb(realScalarObservable * ob,const json_t & j)
{
  return new realScalarEstimator(ob,j);
}

realVectorEstimator* createEstimatorFromOb(realVectorObservable * ob,const json_t & j)
{
  return new realVectorEstimator(ob,j);
}


template<>
struct observableTraits<realScalarObservable >
{
  using storer_t = realScalarStorer;
  using estimator_t = realScalarEstimator;
  
};


template<>
struct observableTraits<realHistogramObservable >
{
  using storer_t = realHistogramStorer;
  using estimator_t = realHistogramEstimator;
  
};


std::vector<estimatorBase*> estimatorFactory:: create(const json_t & j)
  {
    std::vector<estimatorBase*> estimators;
    
    for (const auto & estJson : j )
      {

	std::string id = estJson["kind"].get<std::string>();
	
	if (   estJson.find("forwardWalkingSteps") != estJson.end() )
	  {
	    int i=0;
	    
	    for ( const auto &   fwStepJ : estJson["forwardWalkingSteps"]  )
	      {
		json_t jFW = estJson;
		int steps = fwStepJ.get<int>();
		
		jFW["forwardWalkingSteps"]=steps;
		jFW["targetLabel"]=jFW["label"];
		jFW["label"]=jFW["label"].get<std::string>() + "_fw" + std::to_string(steps);

		if (knownObservableTypes.at(id) == "scalar")
		  {
		    estimators.push_back(new realScalarForwardWalkingEstimator(jFW) );
		  }
		else if (knownObservableTypes.at(id) == "histogram")
		  {
		    estimators.push_back(new realHistogramForwardWalkingEstimator(jFW) );
		  }
		else
		  {
		    throw missingImplementation("Forward walking estimator for " + knownObservableTypes.at(id) + " not yet implemented");
		    
		  }
		i++;
	      }

	   
	    
	  }
	else if (id == "superfluidFraction" )
	  {
	    json_t j2(estJson);
	    
	    auto cmStorer = createSuperfluidFractionStorer(j);
	    const auto & storedSets = cmStorer->sets();
	    std::vector<int> sfSets = j2["sets"].get<std::vector<int> >();
	    
	    if (sfSets.size() != 2)
	      {
		throw invalidInput("Superfluid fraction should have two sets");
	      }
	    int index0,index1;
	    {
	    
	      auto it0 = std::find(storedSets.begin() ,storedSets.end() , sfSets[0] );
	      index0=std::distance(storedSets.begin(), it0);

	      auto it1 = std::find(storedSets.begin() ,storedSets.end() , sfSets[1] );
	      index1=std::distance(storedSets.begin(), it1);
	      
	    }

	    if ( (index0 > storedSets.size()) or (index1 > storedSets.size()) )
	      {
		throw invalidInput("Superfluid states not stored. This is a bug.");
	      }

	    std::vector<int> indices = {index0,index1};
	    j2["sets"]=indices;
	    
	    
	    j2["targetLabel"]="superfluidFractionObservables";
	    
	    estimators.push_back(new superfluidFractionEstimator(j2) );
	    
	    
	  }
	else
	  {
	    
	
	    if ( (id != "forceEnergy") and (id !="energy") )
	      {
		
		estimators.push_back( abstractFactory_t::create(id,estJson) );
		
		
	      }
	  }
      }
    
    return estimators;
  }
