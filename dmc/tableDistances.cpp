#include "tableDistances.h"
#include "geometry.h"
#include "qmcExceptions.h"

void tableDistances::add(int setA)
{
	auto index = _differences.size();
	_differences.resize(_differences.size() + 1);
	_distances.resize(_distances.size() + 1);

	indices1b[setA] = index;
};

void tableDistances::add(int setA,int setB)
{
  if (!isAdded(setA,setB) )
    {
      auto index = _differences.size();
      _differences.resize(_differences.size() + 1);
      _distances.resize(_distances.size() + 1);
  
      indices2b[ std::make_pair(setA,setB)] = index;
    }
};

void tableDistances::add(const std::vector<int> & sets)
{
  if (sets.size() == 1)
    {
      add(sets[0]);
    }
  else if (sets.size() == 2)
    {
      add(sets[0],sets[1]);
    }
  else if (sets.size() ==3)
    {
      add(sets[0],sets[2]);
      add(sets[1],sets[2]);
      add(sets[0],sets[1]);     
      
    }
  
  
};



bool tableDistances::isAdded(int setA,int setB)
{
  return indices2b.find( std::make_pair(setA,setB)  )!=indices2b.end();
};



const tableDistances::diff_t & tableDistances::differences(int setA) const 
{
	int index= indices1b.at(setA);
	return _differences[index];
};

const tableDistances::distances_t & tableDistances::distances(int setA) const 
{
	int index= indices1b.at(setA);
	return _distances[index];
};

const tableDistances::diff_t & tableDistances::differences(int setA, int setB) const 
{
  int index= indices2b.at(std::make_pair(setA,setB));
  return _differences[index];
};

const tableDistances::distances_t & tableDistances::distances(int setA,int setB) const 
{
	int index= indices2b.at(std::make_pair(setA,setB));

	return _distances[index];
};

void tableDistances::update(const tableDistances::states_t & states)
{
	
	for ( const auto  element : indices1b ) // updates one body distances
	{
		_differences[element.second]=geo->differences(states[element.first],{0.,0.,0.});
		_distances[element.second]=norm(_differences[element.second]);

	}

	for ( const auto  element : indices2b ) // updates one body distances
	{
		int setA=element.first.first;
		int setB= element.first.second;
		if (setA == setB)
			{
				_differences[element.second]=geo->differencesTwoBody(states[setA]);
			}
			else
			{
				_differences[element.second]=geo->differencesTwoBody(states[setA],states[setB]);

			}
		_distances[element.second]=norm(_differences[element.second]);

	}	
	
}
