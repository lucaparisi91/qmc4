#include "walkers.h"
#include "energy.h"
#include "wavefunction/productWavefunction.h"
#include "tools.h"
#include <sys/stat.h>

walker::walker(const walker & w)
{
  (*this)=w;
}

void walker::operator=(const walker & w)
{
  _states=w._states; // a vector of particle data
  _tab=w._tab;
  _waveValue=w._waveValue;
  _gradients=w._gradients;
  _lapLog=w._lapLog;
  
  _slaters=w._slaters; // contains the matrix of slater determinants
  
  _phaseGradients=w._phaseGradients;
  storageScalarObservables=w.storageScalarObservables;
  correlatorCurrentTimeIndex=w.correlatorCurrentTimeIndex;
  isFilling=w.isFilling;
}


bool updateForceGradientLaplacian(walker & w,productWavefunction & psi)
{
	/* Update forces ,laplacian and wavefunction value*/
  w.getTableDistances().update(w.getStates());
  w.getTableSlaters().update(w.getStates());
  
  if ( psi.satisfyConstraints(w) )
    {      
      psi.evaluateDerivatives(w);
      return true;
    }
  else
    {
      return false;
    }
};

bool updateForceGradientEnergy(dmcWalker & w,productWavefunction & psi, energy & energyOb)
{
  w.getTableSlaters().update(w.getStates());
  w.getTableDistances().update(w.getStates());
  
  if ( psi.satisfyConstraints(w) )
    {      
      w.getEnergy()=energyOb(w,psi);
      return true;
    }
  else
    {
      return false;
    }
  
};


void dmcWalker::updateMPIDataType()
{
  MPI_Type_free(&getMPIDatatype());
  createMPIDataType();
}

dmcWalker::dmcWalker()
{
  createMPIDataType();
}

dmcWalker::dmcWalker(const dmcWalker & w) : walker::walker(w), _e(w._e)
{
  createMPIDataType();
}

void dmcWalker::operator=(const dmcWalker & w)
{
  walker::operator=(w);
  _e=w._e;
}




void dmcWalker::createMPIDataType()
{
  /*
    Create a data type which send particle configurations, force gradients , laplacian and energy. Does NOT send slater and distances information ( which in standarard dmc are not required for the evolution)
*/
  const int N = getStates().size();
  const int M= 2;
  const int NScalars= 3;
  const int nIntegers = getTimeIndex().size();
  const int nBools = getFillingStatus().size();


  const int NstorageScalarObservables=getStorageScalarCorrelators().size();
  const int Ntot = N*M + NScalars + NstorageScalarObservables + nIntegers + nBools;

  MPI_Aint  offsets[Ntot] ;
  int blockCounts[Ntot];
  MPI_Datatype dtypes[Ntot];
  int k=0;
  for (int i=0;i<N;i++)
    {
      MPI_Get_address(getStates()[i].data(), &offsets[k]);
      blockCounts[k]=getStates()[i].size();
      dtypes[k]=MPI_DOUBLE;
      k++;
    }
  
  for (int i=0;i<N;i++)
    {
      MPI_Get_address(getGradients()[i].data(), &offsets[k]);
      blockCounts[k]=getStates()[i].size();
      dtypes[k]=MPI_DOUBLE;
      k++;
    }
  
  MPI_Get_address(&getLogWave(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;
  k++;
  
  MPI_Get_address(&getLaplacianLog(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;
  k++;
  
  MPI_Get_address(&getEnergy(), &offsets[k]);
  blockCounts[k]=1;
  dtypes[k]=MPI_DOUBLE;  
  k++;
  
  {
  auto it = getStorageScalarCorrelators().begin();
  
  while (it != getStorageScalarCorrelators().end() )
    {
      MPI_Get_address(it->second.data(), &offsets[k]);
      blockCounts[k]=it->second.size();
      dtypes[k]=MPI_DOUBLE;
      k++;
      it++;
    }

  }

  {
  auto it = getTimeIndex().begin();

  while (it != getTimeIndex().end() )
  {
    MPI_Get_address(& (it->second), &offsets[k]);
    blockCounts[k]=1;
    dtypes[k]=MPI_INT;
    k++;
    it++;
  }

  }

  {
  auto it = getFillingStatus().begin();

  while (it != getFillingStatus().end() )
  {
    MPI_Get_address(& (it->second), &offsets[k]);
    blockCounts[k]=1;
    dtypes[k]=MPI_CXX_BOOL;
    k++;
    it++;
  }

  }
    


  
  MPI_Type_create_struct(Ntot, blockCounts, offsets, dtypes, &getMPIDatatype() );
  MPI_Type_commit(&getMPIDatatype());
  
  //MPI_Get_address(getStates().data(), &offsets[0]);
  // MPI_Get_address(&getLogWave() , &offsets[1] );
  // MPI_Get_address(&getGradients() , &offsets[2] );
  // MPI_Get_address(&getLaplacianLog() , &offsets[3] );
  // MPI_Get_address(&getEnergy() , &offsets[4] );  
}




template<class T>
void walkerContainer<T>::push_back(const T &  w)
    {
      _size=_size +1;
    
      if (_size > capacity() )
	{
	  walkers.push_back(std::unique_ptr<T>() );
	  (*(walkers.end() -1 )).reset( new T(w));
	}
      else
	{
	  *(walkers[_size-1])=w;
	}
    }



template<class T>
void walkerContainer<T>::dump(int i)
  {
    std::ofstream f;
    int pId =pTools::rank();

    struct stat st = {0};

    if (stat(baseDir.c_str(), &st) == -1) {
      mkdir(baseDir.c_str(), 0700);
}
    
    f.open(baseDir + "/walkers-Rank" + std::to_string(pId) + ".json");

    json_t j = toJson();
    f << j;
    f.close();
    
  }

template<class T>
json_t walkerContainer<T>::toJson()
  {
    json_t j;
    
    std::vector<json_t > configurations;
    
    for (int i=0;i<walkers.size();i++)
      {
	configurations.push_back(::toJson((*this)[i].getStates()) );
      }
    
    j["configurations"]=configurations;
    
    return j;
  }

template<class T>
void walkerContainer<T>::reserve(size_t size2,const T & w)
  {
    auto oldCap = capacity();
    if (size2 > capacity() ) walkers.resize(size2);
    for (int i=oldCap ; i < capacity() ;i++)
      {
	walkers[i].reset(new T(w));
      }
  }


template<class T>
void walkerContainer<T>::resize(size_t size2)
  {
    if (size2 > capacity()  )
      {
	
	resize(size2,T());
      }
    else
      {
	_size=size2;
      }
  }

template<class T>
void walkerContainer<T>::resize(size_t size2, const T & w)
  {
    reserve(size2,w);    
    _size=size2;
  }






#include "walkers.hpp"
