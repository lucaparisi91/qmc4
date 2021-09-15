#ifndef ABSTRACT_FACTORY_H
#define ABSTRACT_FACTORY_H

#include <map>
#include <iostream>
#include <stdlib.h>
#include <string>
#include "qmcExceptions.h"
#include <string>
#include "traits.h"
using namespace std;

class defaultAction{};

template<class T1,class T2,class creatorType>
class abstractFactory
{
 public:
  typedef T1 object;
  typedef T2 idType;
  
  typedef map<idType,creatorType> creatorMap;
  
  template<class T=defaultAction>
  object* create(const idType &id,const json_t & j , const geometry_t & geo )
  {
    typename creatorMap::const_iterator i;
    i=creators.find(id);
    if (i!= creators.end())
      {
	return (i->second)( j , geo  );
      }
    else
      {
	throw factoryIdNotRecorded(id);
      }
  };



    template<class T=defaultAction>
  object* create(const idType &id,const json_t & j )
  {
    typename creatorMap::const_iterator i;
    i=creators.find(id);
    if (i!= creators.end())
      {
	return (i->second)( j  );
      }
    else
      {
	throw factoryIdNotRecorded(id);
      }
  };

  
  
  bool registerType(const idType &id, creatorType creator )
  {
    return creators.insert(typename creatorMap::value_type(id,creator) ).second;
  };
  
 private:
  
  creatorMap creators;
};

std::string createId(const json_t & j);

#endif
