#include "tools.h"
#include "qmcExceptions.h"
#include  <filesystem>
#include <iostream>
#include <regex>
#include "ptools.h"

std::string ansiColor(const std::string & color_name)
{
	if (color_name == "red") return  "\033[0;31m"; 
	if (color_name == "default") return  "\033[0m";
	if (color_name == "cyan") return  "\033[0;36m"; 
	if (color_name == "green") return  "\033[0;32m"; 
	if (color_name == "green") return  "\033[0;32m"; 
	if (color_name == "yellow") return "\u001b[33m";
	
	throw invalidInput(color_name);

};

bool is_empty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}


json_t toJson(const states_t & states)
{
  json_t jO;
  
  for (int j=0;j<states.size();j++)
    {
      json_t jSet;
      
      for (int i=0;i<getDimensions();i++)
	{
	  const auto & x = states[j].col(i);

	  std::vector<real_t> xCopy(x.data() , x.data() + x.size()) ;
	  
	  if (i==0)
	    jSet["x"]=xCopy;
	  if (i==1)
	    jSet["y"]=xCopy;
	  if (i==2)
	    jSet["z"]=xCopy;
      
      
	}
      jO.push_back(jSet);
    }
  return jO;
}

std::vector<states_t> readStates(const json_t & jI)
{
  std::vector<states_t> states;
  states.resize(jI.size());

  int iStates=0;
  for ( auto & state : states)
    {
      json_t jState=jI[iStates];
      int M = jState.size();
      state.resize(jState.size());
      
      for (int j=0;j<M;j++)
	{
	  json_t jSet = jState[j];

	  int N=jSet["x"].size();
	  auto & confSet = state[j] ;
	  
	  confSet.resize(N,getDimensions() );
	  for (int i=0;i<getDimensions();i++)
	    {
	  
	      std::vector<real_t> x ;
	      
	      if (i==0)
		  x=jSet["x"].get<std::vector<double> >();
	      if (i==1)
		x=jSet["y"].get<std::vector<double> >();
	      if (i==2)
		x=jSet["z"].get<std::vector<double> >();

	      
	      for (int k=0;k<x.size();k++)
		confSet(k,i)=x[k];
      
	    }
	}
      iStates++;
    }
  return states;
}

std::vector<states_t> readStatesFromDirectory(std::string dirname)
{
  //using filesystem = std::filesystem;
  std::vector<states_t> states;
  json_t allStatesJ;
  if (!std::filesystem::exists(dirname) )
    {
     
      return states;
    }

  for (auto & p : std::filesystem::directory_iterator(dirname) )
     {
       std::string filename=p.path();
       std::regex rgx(".*-Rank([0-9]+).*");
       std::smatch match;
       
       bool isMatched=std::regex_match(filename, match, rgx);

       if (isMatched)
	 {
	   json_t jI;
	   std::ifstream f;
	   f.open(filename);
	   f >> jI;
	   allStatesJ.push_back(jI);
	 
	 }
  
	   
     }
  // repeat the number of initial states if smaller than the number of processors

  int i=0;
  int old_size=allStatesJ.size();
  while (allStatesJ.size() < pTools::nProcesses() )
    {
      allStatesJ.push_back(allStatesJ[i%old_size]);
      i++;
    }
  
  // selects the states to load on this process
  int nLocalStates=allStatesJ.size()/pTools::nProcesses();
  
  for (int i=nLocalStates*pTools::rank() ; i< (nLocalStates*(pTools::rank() + 1) );++i)
    {
      const auto & jI=allStatesJ[i];
      auto localStates = readStates(jI["configurations"]);
      states.insert( states.end(), localStates.begin(), localStates.end() );
    };
  
  if (pTools::rank() < (allStatesJ.size() % pTools::nProcesses() ) )
    {
      int i=nLocalStates * pTools::nProcesses() + pTools::rank();
      const auto & jI=allStatesJ[i];
      auto localStates = readStates(jI["configurations"]);
      states.insert( states.end(), localStates.begin(), localStates.end() );
    }
       
       
  
  
  return states;
}

int wrapIndex(int i , int size)
{
  return (i>= 0 ) ? i%size : size  + i;
}



