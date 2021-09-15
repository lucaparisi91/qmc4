#include "abstractFactory.h"

std::string createId(const json_t & j)
  {
    std::string id="";
    if ( j.find("kind") != j.end() )
      id +=j["kind"].get<std::string>();
    
    for (const auto  & element : j)
      {
	if ( element.is_structured() )
	  {
	    auto id2 = createId(element);
	    if (id2 != "")
	      id+="/"+createId(element);
	  }

      }
    return id;
  }

