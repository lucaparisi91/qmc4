#ifndef TABLE_DISTANCES_H
#define TABLE_DISTANCES_H

#include "traits.h"
#include <map>
#include <unordered_map>


class tableDistances
{
public:
  using diff_t = ::difference_t;
  using distances_t = ::distance_t;
  using state_t = ::state_t;
  using states_t = ::states_t;

  tableDistances(){};
  const auto & getGeometry() const {return *geo;}

	void setGeometry(const geometry_t & geo_)  {geo=&geo_;}
	tableDistances(geometry_t & geo_) : geo(&geo_){}

	void update(const states_t & states); // update all registered distances

        void add(int setA );  
        void add(int setA, int setB);

        void add(const std::vector<int> & sets);
  
        bool isAdded(int setA,int setB);
	const diff_t & differences(int setA) const ;
	const diff_t & differences(int setA,int setB) const;
	const auto & differences()  const {return _differences;} 

	const distances_t & distances(int setA) const;
    const distances_t & distances(int setA,int setB) const;

private:
	std::vector<diff_t> _differences; // vectorial distances
	std::vector<distances_t> _distances; // scalar distances
	std::map<std::pair<int,int> , int> indices2b;
	std::unordered_map<int,int> indices1b;
	const geometry_t *  geo;
};

#endif
