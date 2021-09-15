#include "input.h"
#include <regex>




void resolveLink( input && j)
{
  if (! j.is_structured())
    {
      std::string description= j.dump();
      
      std::regex linkExpr("\"@(.*)\"");
      std::smatch matches;

      if ( std::regex_match(description,matches,linkExpr) )
	{
	  std::string filename= matches[1];
	  input newJ=loadFromFile(filename);
	  j=newJ;
	  
	}
    }
  else
    {
      for (auto it= j.begin() ; it!= j.end() ;it++)
	{
	  resolveLink(*it);
	}
    }

  
  
}


void resolveLink( input & j)
{
  resolveLink(std::move(j) );
}





input loadFromFile(std::string filename)
{
  std::ifstream f;
  input j;
  f.open(filename);
  f >> j;
  f.close();
  return j;
};

std::string to_string(const std::string & s)
{
  return s;
  
}


input inputIterator::operator*()
{
  return input(& (*it) );
}

std::istream &  operator>>(std::istream & is , input & inputO)
{  
  is >> (*inputO.jptr() ) ;
  resolveLink(inputO);
  
  return is;
  
};

std::ostream &  operator<<(std::ostream & out , const input & inputO)
{
  out << (*inputO.jptr() ) ;
  return out;
  
}
