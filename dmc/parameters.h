#ifndef PARAMETER_H
#define PARAMETER_H

#include <string>
#include "traits.h"


class optimizationParameter
{
public:
  optimizationParameter(std::string name,size_t length,size_t offset);
  
  size_t leftOffset(){return _offset;}

  size_t size() {return _length; }
  
  
  const std::string & name() const {return _name;}
  
  size_t size() const {return _length;}

private:
  std::string _name;
  size_t _length;
  size_t _offset;
};


class mappedOptimizationParameter : public optimizationParameter
{
public:
  mappedOptimizationParameter(std::string name,size_t length,size_t offset,size_t storageOffset) ;
  
  mappedOptimizationParameter(optimizationParameter parameter,size_t offset);
  
  virtual void setStorageOffset(size_t offset){_storage_offset=offset;};
  
  size_t storageOffset() const {return _storage_offset;}
  
  template<class storage_t>
  auto begin(storage_t & storage) const {return storage.begin() + _storage_offset;}
  
  template<class storage_t>
  auto end(storage_t & storage) const {return storage.begin() + _storage_offset;}
  
  
private:
  size_t _storage_offset;
};

/*
  Contains the unique parameters info. 
  Performs the mapping between a non unique wavefunction mapped parameter
*/
class optimizationParameters
{
public:
  optimizationParameters(){};
  
  bool addParameter(std::string name,size_t size);
  
  const optimizationParameter &  operator[](std::string label);

  size_t size() const {return _parameters.size() ;}
  
  mappedOptimizationParameter mapParameter(optimizationParameter & param,std::string label);
  
private:
  int currentOffset=0;
  std::vector<mappedOptimizationParameter> _parameters;
  std::map<std::string,int > nameToIndexMap;
};

class wavefunctionParameters
{
  public:
    wavefunctionParameters(){};

    void addParameter(int wave,optimizationParameter parameter, std::string key); // register an optimization parameter


    const auto & parameters() const {return _parameters;}

    size_t size() const {return parametersInfo.size();}


  private:
    std::vector<std::vector<mappedOptimizationParameter>  > _parameters;

    optimizationParameters parametersInfo; 
};


#endif
