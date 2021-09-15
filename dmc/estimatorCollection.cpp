#include "estimatorCollection.h"
#include "estimators.h"
#include "correlationEstimator.h"


void estimatorCollection::clear()
	{
		for (auto & est : _estimators)
		{
			est->clear();
		}
	}


void estimatorCollection::dump()
		{
			for (auto & est : _estimators)
			{
				est->dump();
			}
		}


void estimatorCollection::accumulateMPI(int root)
		{
			for (auto & est : _estimators)
			{
			  est->accumulateMPI(root);
			}
		}

void estimatorCollection::accumulate( estimatorCollection::walker_t & w,estimatorCollection::wavefunction_t & psi)
{
  
  for (auto & storerPtr : _storers)
    {
      storerPtr->store(w,psi);
    }
  for (auto & est : _estimators)
    {
      est->accumulate(w,psi);
    }
}


void estimatorCollection::reserve( walker_t & w)
{
  for (auto & storerPtr : _storers)
    {
      storerPtr->reserve(w);
    }
}

