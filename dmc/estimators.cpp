#include "estimators.h"
#include <fstream>
#include "tools.h"
#include "ptools.h"
#include <iomanip>


template<class observable_t>
void estimator<observable_t>::write(std::ostream & stream)
{
  acc.writeAverage(stream) ;
};


template class estimatorObservable<realScalarObservable>;

estimatorBase::estimatorBase(std::string label) : _label(label)
{
	filename= label+".dat";
	
	f.open(filename,std::fstream::out | std::fstream::app);
}

estimatorBase::~estimatorBase()
{
	f.close();
}

void estimatorBase::dump()
{
	write(f);
	f << std::endl;
}

realScalarEstimator::realScalarEstimator(std::string label_,realScalarObservable * ob_) : estimatorObservable<realScalarObservable>::estimatorObservable(label_,ob_)
{
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << getLabel() << std::endl;
    }
  of.close();
}

realScalarEstimator::realScalarEstimator(realScalarObservable * ob_,const json_t & j) : realScalarEstimator(j["label"].get<std::string>(),ob_ ){}



realHistogramEstimator::realHistogramEstimator(std::string label,realHistogramObservable * ob_,size_t size,real_t minx,real_t maxx) : estimatorObservable<realHistogramObservable>::estimatorObservable(label,ob_)
{
  getAccumulator()=accumulator_t(size,minx,maxx);
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x " << getLabel() << std::endl;
    }
  
  of.close();

  x.resize(size);
  for (int i=0;i<x.size();i++)
    {
      x[i]=getAccumulator().stepSize() * i;
    }
  
}

realHistogramEstimator::realHistogramEstimator(realHistogramObservable * ob_,const json_t & j): realHistogramEstimator(j["label"].get<std::string>(),ob_,j["bins"].get<size_t>(), j["minx"].get<real_t>(),j["maxx"].get<real_t>()   )
{
  
}


void realHistogramEstimator::write(std::ostream & stream)
{
  auto out=getAccumulator().average();
  for (int i=0;i<x.size();i++)
    {
      stream << x[i] << " " <<  out[i] << std::endl ;
    }
};

void realScalarEstimator::write(std::ostream & stream)
{
    stream << std::scientific <<  std::setprecision(12);

  if ( getAccumulator().getWeight() > 0 )
    {
      estimatorObservable<realScalarObservable>::write(stream);
    }
  
}


void realVectorEstimator::write(std::ostream & stream)
{
  auto out=getAccumulator().average();
  auto & x = getObservable().x();
  stream << std::scientific <<  std::setprecision(12);

  for (int i=0;i<x.size();i++)
    {
      stream << x[i] << " " <<  out[i] << std::endl ;
    }
};



realVectorEstimator::realVectorEstimator(std::string label,realVectorObservable * ob_) : estimatorObservable<realVectorObservable>::estimatorObservable(label,ob_)
{
  auto size= getObservable().x().size();
  
  getAccumulator()=accumulator_t(size);
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x " << getLabel() << std::endl;
    }
  
  of.close();
  
}

realVectorEstimator::realVectorEstimator(realVectorObservable * ob_,const json_t & j) : realVectorEstimator::realVectorEstimator(j["label"].get<std::string>(),ob_)
{
  
}

