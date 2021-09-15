#include "correlationEstimator.h"
#include "walkers.h"
#include "tools.h"

realScalarStorer::realScalarStorer(std::string label_,realScalarObservable * ob_, int recordSteps_) : ob(ob_),recordSteps(recordSteps_),storer::storer(label_) {}


realScalarStorer::realScalarStorer(realScalarObservable * ob_,const json_t & j )  : realScalarStorer(j["label"].get<std::string>(),ob_,j["recordSteps"].get<int>()  ) {};

void realScalarStorer::reserve(realScalarStorer::walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps);
  w.getTimeIndex()[getLabel()]=0.;
  w.getFillingStatus()[getLabel()]=true;
}

void realScalarStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
  w.getFillingStatus().at(getLabel())=true;
}


void realScalarStorer::store( walker_t & w, wavefunction_t & psi )
  {
    int & i = w.getTimeIndex().at(getLabel() );
    
    w.getStorageScalarCorrelators().at(getLabel() )(i)=(*ob)(w,psi);
    i=(i+1)% recordSteps;
    if (i==0)
      {
	w.getFillingStatus().at(getLabel())=false;
      }
  }

void realScalarForwardWalkingEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    const auto & i = w.getTimeIndex().at(targetLabel);
    
    auto & data =w.getStorageScalarCorrelators().at(targetLabel);
    auto recordSteps = data.size();
    assert(forwardWalkingSteps < recordSteps);
    auto j =  wrapIndex(i - forwardWalkingSteps - 1 ,recordSteps ) ;
    const auto & v = data( j);
    
    if ( ! w.getFillingStatus().at(targetLabel) )
      {
	getAccumulator()+=v;
      }
    
  }

realScalarForwardWalkingEstimator::realScalarForwardWalkingEstimator(std::string label,std::string targetLabel_, int forwardWalkingSteps_) : estimator<realScalarAccumulator_t>(label),targetLabel(targetLabel_),forwardWalkingSteps(forwardWalkingSteps_)
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

realScalarForwardWalkingEstimator::realScalarForwardWalkingEstimator( const json_t & j) : realScalarForwardWalkingEstimator( j["label"].get<std::string>() ,  j["targetLabel"].get<std::string>(), j["forwardWalkingSteps"].get<int>() )
{
  
}

void realScalarForwardWalkingEstimator::write(std::ostream & stream)
{
    if ( getAccumulator().getWeight() > 0 )
    {
      estimator<realScalarAccumulator_t>::write(stream);
    }

}



realHistogramStorer::realHistogramStorer(std::string label_,realHistogramObservable * ob_,size_t size,real_t minx,real_t maxx,int recordSteps_) : ob(ob_),recordSteps(recordSteps_),storer::storer(label_) ,tmpAcc(size,minx,maxx)
{
  
}

void realHistogramStorer::store(realHistogramStorer::walker_t & w,wavefunction_t & psi)
{
  int & i = w.getTimeIndex().at(getLabel() );
  tmpAcc.clear();
  ob->accumulate(w,psi,tmpAcc);
  
  auto & data=w.getStorageScalarCorrelators().at(getLabel() );
  
  for(int j=0;j<tmpAcc.size();j++)
    {
      data(i *tmpAcc.size() + j)= tmpAcc.sums()(j)/tmpAcc.weight();
    }
  
  i=(i+1)% recordSteps;
  if (i==0)
    {
      w.getFillingStatus().at(getLabel())=false;
    }
}


void realHistogramStorer::reserve(walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps*tmpAcc.size() );
  w.getTimeIndex()[getLabel()]=0.;
  w.getFillingStatus()[getLabel()]=true;
  
}

void realHistogramStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
  w.getFillingStatus().at(getLabel())=true;
}

realHistogramStorer::realHistogramStorer(realHistogramObservable * ob_,const json_t & j ) : realHistogramStorer(j["label"].get<std::string>(),ob_,j["bins"].get<int>(),j["minx"].get<real_t>(),j["maxx"].get<real_t>(),j["recordSteps"].get<int>())
{
  
}

void realHistogramForwardWalkingEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    const auto & i = w.getTimeIndex().at(targetLabel);

    auto & acc=getAccumulator();
    
    
    auto & data =w.getStorageScalarCorrelators().at(targetLabel);
    auto recordSteps = data.size()/getAccumulator().size();
    assert(forwardWalkingSteps*acc.size() < data.size());
    
    auto j =  wrapIndex(i - forwardWalkingSteps - 1 ,recordSteps ) ;
    
    auto & sums = acc.sums();
    if ( ! w.getFillingStatus().at(targetLabel) )
      {
	for(int ii=0;ii<acc.size();ii++)
	  {
	    sums(ii)+=data(j*acc.size() + ii );
	  }
	acc.weight()+=1;
      }
    
  }



realHistogramForwardWalkingEstimator::realHistogramForwardWalkingEstimator(std::string label,std::string targetLabel_,size_t size,real_t minx,real_t maxx,int forwardWalkingSteps_) :  estimator<realHistogramAccumulator_t>(label),targetLabel(targetLabel_),forwardWalkingSteps(forwardWalkingSteps_)
{
  getAccumulator().resize(size,minx,maxx);


  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x" <<" "<<getLabel() << std::endl;
    }
  
  of.close();
}


void realHistogramForwardWalkingEstimator::write(std::ostream & stream)
{
    if ( getAccumulator().weight() > 0 )
    {
      estimator<accumulator_t>::write(stream);
    }

}

realHistogramForwardWalkingEstimator::realHistogramForwardWalkingEstimator(const json_t & j) : realHistogramForwardWalkingEstimator(j["label"].get<std::string>(),j["targetLabel"].get<std::string>(),j["bins"].get<int>(),j["minx"].get<real_t>(),j["maxx"].get<real_t>(),j["forwardWalkingSteps"].get<int>() )
{
  
}


centerOfMassStorer::centerOfMassStorer(std::string label_,int recordSteps_,std::vector<int> sets_) :
  storer::storer(label_),
  recordSteps(recordSteps_),
  _sets(sets_)
{
  xcm.resize(_sets.size()*getDimensions() );
  
}


centerOfMassStorer::centerOfMassStorer(const json_t & j) :
  centerOfMassStorer::centerOfMassStorer(
					 j["label"].get<std::string>(),
					 j["bins"].get<int>(),
					 j["sets"].get<std::vector<int> >()
)
{
  
}

void centerOfMassStorer::reserve(walker_t & w)
{
  w.getStorageScalarCorrelators()[getLabel()].resize(recordSteps*xcm.size() );
  w.getStorageScalarCorrelators()[getLabel()].setConstant(0);
  w.getTimeIndex()[getLabel()]=0.;
  w.getFillingStatus()[getLabel()]=true;
}

void centerOfMassStorer::computeCenterOfMass(walker_t & w)
{
  std::fill(xcm.begin(),xcm.end(),0);
  int k=0;

  
  for ( auto set : _sets)
    {
      auto & positions = w.getStates()[set];
      
      const int N = getN(positions);
      
      for (int i=0;i<N;i++)
	{
	  for (int d=0;d<getDimensions();d++)
	    {
	      xcm[k+d]+=positions(i,d);
	    }
	  
	}
      ;
      k+=getDimensions();
    }

  
}


void centerOfMassStorer::store(realHistogramStorer::walker_t & w,wavefunction_t & psi)
{
  int & i = w.getTimeIndex().at(getLabel() );
  
  computeCenterOfMass(w);
  
  auto & data=w.getStorageScalarCorrelators().at(getLabel() );
  
  for(int j=0;j<xcm.size();j++)
    {
      data(i *xcm.size() + j)= xcm[j];
    }
  
  i=(i+1)% recordSteps;
  if (i==0)
    {
      w.getFillingStatus().at(getLabel())=false;
    }
  
}

void centerOfMassStorer::reset(walker_t & w)
{
  w.getTimeIndex().at(getLabel() )=0.;
  w.getFillingStatus().at(getLabel())=true;
}


superfluidFractionEstimator::superfluidFractionEstimator(std::string label,std::string targetLabel_, size_t size,int setA_, int setB_) :
  estimator<realVectorAccumulator_t>(label),
  targetLabel(targetLabel_),
  setA(setA_),
  setB(setB_)
{
  getAccumulator().resize(size);
  
  std::ifstream of;
  of.open(getFileName());
  if (is_empty(of) & pTools::rank() == 0 )
    {
      auto & f=getFileDescriptor();
      f << "x" <<" "<<getLabel() << std::endl;
    }
  
  of.close();
}


superfluidFractionEstimator::superfluidFractionEstimator(const json_t & j) :
  superfluidFractionEstimator::superfluidFractionEstimator( 
							   j["label"].get<std::string>(),
							   j["targetLabel"].get<std::string>(),
							   j["bins"].get<int>(),
							   j["sets"][0].get<int>(),
							   j["sets"][1].get<int>()
							    )
  
{
  
}

void superfluidFractionEstimator::write(std::ostream & stream)
{
  int i=0;
  const auto & sums = getAccumulator().sums();
  const auto & weights = getAccumulator().weights();
  
  for (size_t i=0;i<sums.size();i++)
    {
      if ( weights(i) > 0 )
	{
	  stream << i <<  " " << sums(i)/(weights(i) ) << std::endl;
	}
    }
  
}


void superfluidFractionEstimator::accumulate(walker_t & w,wavefunction_t & psi)
  {
    
    const auto & i = w.getTimeIndex().at(targetLabel);

    
    auto & acc=getAccumulator();
    
    const Eigen::ArrayXd & data=w.getStorageScalarCorrelators().at(targetLabel );
    
    
    size_t recordSize = data.size()/getAccumulator().size();
    
    int j =  wrapIndex(i - 1 ,acc.size() ) ;
    
    int recordOffsetA=getDimensions() * setA;
    int recordOffsetB=getDimensions() * setB;

    real_t W=0;
    //std::cout << j << std::endl;
    //std::cout << recordSize << std::endl;

    //std::cout << data.size() << std::endl;
    for (int i=0;i<=j;i++)
      {
	W=0;
	for (int id=0;id<getDimensions();id++)
	  {
	    real_t x1id=data(i*recordSize + recordOffsetA + id);
	    real_t x1jd=data(j*recordSize + recordOffsetA + id);
	    real_t x2id=data(i*recordSize + recordOffsetB + id);
	    real_t x2jd=data(j*recordSize + recordOffsetB + id);
	    
	    W+=(x1id - x1jd)*(x2id-x2jd);
	    
	  }
	acc.accumulate(W,j-i);
	
      }
    
  }



centerOfMassStorer* createSuperfluidFractionStorer(const json_t & j)
{
  json_t j2;
  bool foundSF=false;
  
  std::vector<int> sets;
   for (const auto & storeJson : j )
      {
	std::string id = storeJson["kind"].get<std::string>();

	if (id == "superfluidFraction")
	  {
	    if (not foundSF )
	      {
		j2=storeJson;
		j2["label"]="superfluidFractionObservables";
		foundSF=true;
	      }
	    
	    for (auto set : storeJson["sets"].get<std::vector<int> >() )
	      {
		if ( std::find(sets.begin(),sets.end(),set) == sets.end() )
		  {
		    sets.push_back(set);
		  }
	      }
	    
	    
	  }
      }

   if (foundSF)
     {
       j2["sets"]=sets;
       return new centerOfMassStorer(j2);
     }
   else
     {
       return NULL;
     }
   
}
