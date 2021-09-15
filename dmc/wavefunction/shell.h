#ifndef SHELL_H
#define SHELL_H

#include <complex>
#include <cmath>
#include <vector>
#include <tuple>
#include <iostream>
#include <algorithm>


#if DIMENSIONS == 3


class shell
{
public:
  shell(unsigned int n_) : n(n_){list();};
  size_t capacity() const {return shellStructure.size();}
  unsigned int level() const {return n;}
  const auto & operator[](unsigned int i) const{return shellStructure[i] ;}
 private:
  unsigned int n;
  void list()
  {
    int nMax=int(sqrt(n));
    
    for (int nx=-nMax;nx<=nMax;nx ++ )
      {
	for (int ny=-nMax;ny<=nMax;ny++ )
	  {
	    for (int nz=-nMax;nz<=nMax;nz++ )
	      {
		if (nx*nx + ny*ny + nz*nz==n)
		  {
		    shellStructure.push_back(std::make_tuple(nx,ny,nz) );
		   
		  }
		
		
	      }

	  }
      }
      
  }
  
  std::vector<std::tuple<int, int,int> > shellStructure;
  
};


class shellStructure
{
public:
  using shell_t = shell;
  
  shellStructure(unsigned int N) : _N(N){
    
    shells.push_back(shell_t(0));
    fillingNumbers.push_back(1);
    int n=1;
    while(fillingNumbers.back() < N)
      {
	auto sh=shell_t(n);
	if ( sh.capacity()>0)
	  {
	    shells.push_back(sh);
	    fillingNumbers.push_back(shells.back().capacity());
	    fillingNumbers.back()+=fillingNumbers[fillingNumbers.size()-2];
	    
	  }
	n++;
	
      }
    
  };
  
  bool isMagicNumber(int n)
  {
    return not (std::find(fillingNumbers.begin(),fillingNumbers.end(),n)==fillingNumbers.end());
  }
  
  const auto & getMagicNumbers() const {return fillingNumbers;}

  const shell_t &  getShell(size_t i) const {return shells[i];}

  const shell_t & operator[](size_t i) const {return getShell(i);}

  size_t nShells() const {return shells.size();}

  size_t capacity() const
  {
    return fillingNumbers[-1];
  };
  
private:
  int _N;
  std::vector<shell_t> shells;
  std::vector<int> fillingNumbers;
  
  std::vector<std::tuple<int,int> > orderedShells;
};


std::ostream & operator<<(std::ostream & out,const shell & shell);

std::ostream & operator<<(std::ostream & out,const shellStructure & sh);

#endif


#endif
