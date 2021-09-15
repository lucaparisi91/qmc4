#include "ptools.h"
#include <iostream>
#include <numeric>
#include "walkers.h"

namespace pTools
{

  int init(int argc,char** argv) {return MPI_Init(&argc,&argv);}
  
  int finalize(){return MPI_Finalize();}

  int nProcesses() {int res;  MPI_Comm_size(MPI_COMM_WORLD, &res); return res; }

  int rank() {int res;  MPI_Comm_rank(MPI_COMM_WORLD, &res); return res; }
  
  bool isMaster()
  {
    if (rank()==0)
      return true;
    else
      return false;
  }
  

  int broadcast(int * data, int count,int root)
  {
    return MPI_Bcast( data, count,MPI_INT,root,MPI_COMM_WORLD);
  }
  int broadcast(char * const data, int count,int root)
  {
    return MPI_Bcast( data, count,MPI_CHAR,root,MPI_COMM_WORLD);
  }

  
  int broadcast(int * data,int root)
  {
    return broadcast(data,1,root);
  }

  int broadcast(std::string * data,int root)
  {
    int size=data->size();
    broadcast(&size,root);
    data->resize(size);
    return broadcast(const_cast<char*>(data->data()),size,root); 
  }

  
  int sum(double * send,double * recv ,int count,int root)
  {
    return MPI_Reduce( send, recv,count, MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
  };

  double sum(const double & sum, int root )
  {
    double tmp=sum;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }

  double sumAll(const double & sum)
  {
    double tmp;
    int status = MPI_Allreduce( & sum, & tmp, 1,
			       MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return tmp;
  }

  int sumAll(const int & sum)
  {
    int tmp;
    int status = MPI_Allreduce( & sum, & tmp, 1,
			       MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    return tmp;
  }

  
  int sum(const int & sum, int root )
  {
    int tmp=sum;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_INT,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }
  
  size_t sum(const size_t & sum, int root )
  {
    size_t tmp=0;
    int status = MPI_Reduce (&sum,&tmp,1,MPI_SIZE_T,MPI_SUM,root,MPI_COMM_WORLD);
    return tmp;
  }

  void determineLoadBalanceComunicationsAliasMethod( std::vector<int> & populations,std::vector<int> & permutations,std::vector<int> & sources,std::vector<std::vector<int> > & destinations,std::vector<int> & amounts)
  {
  /* 
Sources: filled with the rank of the processor it must receive from
amounts: filled with number of walkers send by the source
populations : current population distribution
  */
  
  int size = std::accumulate(populations.begin(),populations.end(),0);
  permutations.resize(size);
  std::iota (std::begin(permutations), std::end(permutations), 0);
  // sort such that A[i] < k all on the left of A[j] > k
  int k= size/populations.size();
  // add phantom tasks
  int r=size%populations.size();
  
  int i=0;
  while( (i < populations.size())  and ( populations[i] <= k) )
    {
      i+=1;
    }
  
  int j=i+1;
  while (j<populations.size() )
    {
      if (populations[j] <= k )
	{
	  std::swap(populations[i],populations[j]);
	  std::swap(permutations[i],permutations[j]);
	  i+=1;
	}
      j+=1;
      
    }
  
  // alias algorithm to determine the amount of comunications
  
  j=0;
  sources.resize(populations.size(),0);
  amounts.resize(populations.size(),0);
  destinations.resize(0);
  destinations.resize(populations.size(),{});
  std::fill(sources.begin(),sources.end(),0);
  std::fill(amounts.begin(),amounts.end(),0);
  
  
  
  int I,J;
  auto balanced =  [&] ( int idx) {return  k + (1 ? idx < r : 0 );  } ;
  
  
  while ( (i> j) and (i< populations.size()) )
    {
      J = permutations[j];
      I = permutations[i];
      
      sources[J]=I;
      
      amounts[J]=balanced(j)- populations[j];
      populations[i]=populations[i] - amounts[J];
      destinations[I].push_back(J);
      
      if (populations[i] < balanced(i) )
	{
	  ++i;
	}
      j++;
    }
  }  


walkerDistribution::walkerDistribution()
{
  _currentRank=rank();
  _nProcesses=nProcesses();
  walker_tag=500;
  populations.resize(_nProcesses);
  nWalkersReceived.resize(_nProcesses);
  _permutations.resize(_nProcesses);
  tmpPopulations.resize(_nProcesses);
  _sources.resize(_nProcesses);
  _sendToRanks.resize(_nProcesses);
  inFlight=false;
}
  
void walkerDistribution::determineComm(const std::vector<int> & populations)
{
  tmpPopulations=populations;
  
  pTools::determineLoadBalanceComunicationsAliasMethod( tmpPopulations, _permutations, _sources, _sendToRanks, nWalkersReceived);

  // if (rank() == 0)
  //   {
      
  //     for (int i=0;i<_sendToRanks.size();i++)
  // 	{
  // 	  for(int j=0;j<_sendToRanks[i].size();j++)
  // 	    {
  // 	      int dest=_sendToRanks[i][j];
  // 	      std::cout << i << " ->" << dest << " " << nWalkersReceived[dest] << std::endl;
  // 	    }
  // 	}

      
  //   }
    }

  int isend(double * p,int count,int dest,int tag,MPI_Request &req)
  {
    return MPI_Isend( p ,  count,MPI_DOUBLE, dest, tag,
			  MPI_COMM_WORLD, &req);
  }

  
  int isend(state_t & state,int toRank,int tag,MPI_Request &req)
  {
    auto data_ptr = state.data();
    return isend(data_ptr,state.size(),toRank,tag,req);
  };

  int partialSend(dmcWalker & w,int destination,int tag)
  {
    w.updateMPIDataType();
    return MPI_Send(MPI_BOTTOM, 1, w.getMPIDatatype(), destination, tag, MPI_COMM_WORLD);
  }
  
  int partialRecv(dmcWalker * w, int source,int tag)
  {
    w->updateMPIDataType();
    MPI_Status status;
    
    return MPI_Recv(MPI_BOTTOM, 1, w->getMPIDatatype(), source, tag, MPI_COMM_WORLD, & status);

  }

  
  int ipartialSend(dmcWalker & w,int destination,int tag, MPI_Request * req)
  {
    w.updateMPIDataType();
    return MPI_Isend(MPI_BOTTOM, 1, w.getMPIDatatype(), destination, tag, MPI_COMM_WORLD,req);
  }
  
  int ipartialRecv(dmcWalker * w, int source,int tag,MPI_Request * req)
  {
    w->updateMPIDataType();
    MPI_Status status;
    
    return MPI_Irecv(MPI_BOTTOM, 1, w->getMPIDatatype(), source, tag, MPI_COMM_WORLD, req);
  }

  
    std::vector<int> walkerDistribution::gatherPopulations(int localPopulation)
  {
    std::vector<int> _populations;
    _populations.resize(_nProcesses);
    // every process contains the current population size
    MPI_Allgather(
		  & localPopulation,
		  1,
		  MPI_INT,
		  _populations.data(),
		  1,
		  MPI_INT,
		  MPI_COMM_WORLD);
    
    // if (rank()==0)
    //   {
    // 	for (auto &p : _populations)
    // 	  {
    // 	    std::cout << p <<" " ;
    // 	  }
    // 	std::cout << std::endl;
    //   }

    return _populations;  
  }


void walkerDistribution::isendReceive(walkerDistribution::walkers_t & walkers)
  {
    // all gather population ditribtuion accross processors
    populations=gatherPopulations(walkers.size() );
    // determine what to send to whom
    tmpPopulations=populations;
    determineComm(tmpPopulations);
    
    // // non blocking send
     auto & destinations = _sendToRanks[_currentRank];
    sendRequests.resize(0);
    localSentWalkers=0;
    int k=0;
    
     for (int i=0;i<destinations.size();i++)
       {
     	int dest = destinations[i];
	
     	const auto & amount = nWalkersReceived[dest];
   	localSentWalkers+=amount;

	 // if (pTools::rank() == 0)
	 //  {
	 //    std::cout << _currentRank  << "-> " << dest << " " << amount << std::endl;
	 //  }
	
     	sendRequests.resize(sendRequests.size() + amount);
     	for(int j=0;j<amount;j++)
     	  {
	    
     	    ipartialSend(walkers[walkers.size() -1 - k  ],dest,walker_tag + j,&sendRequests[k]);
     	    k++;
	  }
	
       }
     
     // non blocking receive
    
     const auto & amount = nWalkersReceived[_currentRank];
     
     if (amount!=0)
       {
	  // if (pTools::rank() == 0 )
	  //   {
	  //     std::cout << _currentRank << "<-" << _sources[_currentRank] << " " << amount << std::endl;
	  //   }
	 
     	 walkers.reserve(walkers.size() + amount , walkers[std::max(walkers.size()-1,(size_t)0)]);
	 receiveRequests.resize(amount);
	 
	 for( int j=0;j < amount ; j++ )
	   {
	     ipartialRecv(& walkers[walkers.size() + j   ],_sources[_currentRank],walker_tag + j ,&receiveRequests[j]);
       }

       }
     walkers.resize(walkers.size() - localSentWalkers);
     inFlight=true;
     
  }

int walkerDistribution::wait(walkerDistribution::walkers_t & walkers)
{
  
  if (inFlight)
    {
     
      auto & amount = nWalkersReceived[_currentRank];
  
  // if( pTools::rank() == 0)
  //   {
  //     std::cout << "lw: " << localSentWalkers << std::endl;
  //     std::cout << "cap: " << walkers.capacity() << std::endl;
  //     std::cout << "amount: " << amount << std::endl;
  //     std::cout << "size: " << walkers.size() << std::endl;
  //     std::cout << "receive_from "  <<  _sources[_currentRank]<<std::endl<<std::flush;
  //     auto & destinations = _sendToRanks[_currentRank];
  //     for (int i=0;i<destinations.size();i++)
  // 	{
	  
  // 	  std::cout << _currentRank << "!->" << destinations[i] << std::endl<<std::flush;
  // 	}
  //   }
  // MPI_Barrier(MPI_COMM_WORLD);
      
      
      
  if (amount != 0)
    {
      waitStatusReceives.resize(receiveRequests.size());
      MPI_Waitall(receiveRequests.size(),receiveRequests.data(), waitStatusReceives.data());
      
    }
  
  if ( sendRequests.size() > 0 )
    {
      waitStatusSends.resize(sendRequests.size() );
      MPI_Waitall( sendRequests.size(),sendRequests.data(),waitStatusSends.data());
    }

  // if (pTools::rank() == 1)
  //   {
  //     std::cout << "localSentWalkers: " << localSentWalkers << std::endl;
  //     std::cout << "nAmount: " << amount << std::endl;    
  //     std::cout << "last: " << walkers.size() + localSentWalkers + amount  << std::endl;    
  //     std::cout << "eLast: " << walkers[ walkers.size() + localSentWalkers + amount - 1].getEnergy() << std::endl;    
  //   }
  
      auto i2= walkers.size() + localSentWalkers;
      auto n = std::min(localSentWalkers,   amount );
      auto it2 = walkers.begin() + i2;
      auto it = walkers.begin() + walkers.size();

      
      for ( int i=0;i<n;i++ )
     {
       //std::cout << (*it2)->getEnergy() << std::endl;
       std::swap(*it2,*it);
       it2++;
       it++;
     }

      
      walkers.resize(walkers.size() + amount);

      inFlight=false;
    }
  
  return 0;

}



int sum(  Eigen::Matrix<real_t,Eigen::Dynamic,1> & vec, int root )
  {
    int err;
    if (rank() == root)
      {
	err =MPI_Reduce (MPI_IN_PLACE,vec.data(),vec.rows(),MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
      }
    else
      {
	err =MPI_Reduce (vec.data(),vec.data(),vec.rows(),MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
      }

    return err;
   
  }

int sum(  Eigen::Array<real_t,Eigen::Dynamic,1> & vec, int root )
  {
   
     int err;
    if (rank() == root)
      {
	err =MPI_Reduce (MPI_IN_PLACE,vec.data(),vec.rows(),MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
      }
    else
      {
	err =MPI_Reduce (vec.data(),vec.data(),vec.rows(),MPI_DOUBLE,MPI_SUM,root,MPI_COMM_WORLD);
      }

    return err;
  }


  
};

