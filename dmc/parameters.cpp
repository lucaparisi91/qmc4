#include "parameters.h"
#include "qmcExceptions.h"

optimizationParameter::optimizationParameter(std::string name,size_t length,size_t offset) : _name(name),_length(length),_offset(offset) {}


bool optimizationParameters::addParameter(std::string name,size_t size)
{
  bool nameExists=std::find_if(_parameters.begin(), _parameters.end() , [&name](optimizationParameter & param) {return param.name() == name;} )!= _parameters.end();
  
  if ( !nameExists)
    {
      mappedOptimizationParameter  param(name,size,0,currentOffset);
      _parameters.push_back(param);
      currentOffset+=size;
      nameToIndexMap[name]=_parameters.size()-1;
      
      return true;
    }
  else
    {
      return false;
    }
  
}


mappedOptimizationParameter::mappedOptimizationParameter(std::string name,size_t length,size_t offset,size_t storageOffset) : optimizationParameter::optimizationParameter(name,length,offset) ,_storage_offset(offset) {}


mappedOptimizationParameter::mappedOptimizationParameter(optimizationParameter parameter, size_t offset) : optimizationParameter::optimizationParameter(parameter) , _storage_offset(offset) {}


mappedOptimizationParameter optimizationParameters::mapParameter(optimizationParameter & param,std::string label)
{
  auto it = nameToIndexMap.find(label);
  
  if ( it != nameToIndexMap.end() )
    {
      mappedOptimizationParameter mappedParam(param, _parameters[it->second].storageOffset() );
      
      return mappedParam;      
    }
  else
    {
      throw  invalidInput("Label not found");
    }
}
