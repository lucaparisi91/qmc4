#ifndef TOOLS_H
#define TOOLS_H

#include "traits.h"
#include <fstream>
#include <vector>
#include <iostream>

json_t toJson(const states_t & states);
std::vector<states_t> readStates(json_t & jI);

std::string ansiColor(const std::string & color_name);

bool is_empty(std::ifstream& pFile);


inline int getN(const state_t & state){  return state.rows();};

inline size_t size(const distance_t & dis) {return dis.rows();}

inline constexpr int getDimensions() { return DIMENSIONS; }

std::vector<states_t> readStatesFromDirectory(std::string filename);



int wrapIndex(int i , int size);

//inline constexpr int getDimensions(const state_t & state) {return DIMENSIONS;}


template<class T>
double findRootBrente(T f,double a,double b,double eps)
{
  
  if (f(a)*f(b)>=0)
    {
      std::cout << "Not alternating ends!";
      exit(0);
      return a;
    }
  
  double c,s,tmp,d;
  bool flag;
  if (std::abs(f(a)) < std::abs(f(b)) )
    {
      // swap a and b
      tmp=a;
      a=b;
      b=tmp;
    }
  
  c=a;
  s=a;
  d=a;
  
  flag=true;
  
  while( ! (f(b)==0 or f(s)==0 or abs(b-a)<=eps))
    {
      if (  (f(a) != f(c)) and (f(b)!=f(c)))
	{
	  // use quadratic interpolation to find the closest root
	  s=a*f(b)*f(c)/( (f(a) - f(b))*(f(a)-f(c)) ) + b*f(a)*f(c)/((f(b)-f(a))*( f(b)-f(c) )) + c*f(a)*f(b)/( (f(c)-f(a))*( f(c)-f(b) ) );
	}
      else
	{
	  // metodo della secante
	  s=b-f(b)*(b-a)/(f(b)-f(a));
	}
      if (
	  (s < (3*a+b)/4 or s>b) or
	  (flag and (abs(s-b) >= abs(b-c)/2)) or
	  (! flag and (abs(s-b)>=abs(c-d)/2)  ) or
	  (  flag and (abs(b-c) < eps  ) ) or
	  ( !flag and (abs(c-d) < eps)  ) 
	  )
	{
	  s=(a+b)/2.;
	  flag=true;
	}
      else
	{
	  flag=false;
	}
      d=c;
      c=b;
      if (f(a)*f(s)<0)
	{
	  b=s;
	}
      else
	{
	  a=s;
	}

      if (abs(f(a))<abs(f(b)))
	{
	  // swap a and b
	  tmp=a;
	  a=b;
	  b=tmp;
	}
	
    }
  
  if (f(s)==0)
    {
      return s;
    }
  else
    {
      return b;
    }
}







#endif
