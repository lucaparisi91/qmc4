#include "shell.h"

#if DIMENSIONS == 3

std::ostream & operator<<(std::ostream & out,const shell & shell)
{
  std::cout << "Shell: "<<shell.level()<< " Capacity: "<<shell.capacity()<<std::endl;

  for(int i=0;i<shell.capacity();i++)
    {
      out << std::get<0>(shell[i]) << " " << std::get<1>(shell[i]) << " " << std::get<2>(shell[i])<<std::endl;
	
    }
  return out;
}


std::ostream & operator<<(std::ostream & out,const shellStructure & sh)
{
  out << "Magic numbers: ";
  
  for(int i=0;i<sh.getMagicNumbers().size();i++)
    {
      out << sh.getMagicNumbers()[i] << ",";
    }
  out << std::endl;
  
  for (int i=0;i<sh.nShells();i++)
    {
      out << sh.getShell(i) << std::endl;
    }
  
  return out;
}



#endif
