#include "pimcConfigurations.h"
#include "qmcExceptions.h"
#include <filesystem>
#include <sstream>
#include <algorithm>
#include "hdf5IO.h"
#include <iostream>

namespace fs = std::filesystem;

namespace pimc
{
    
int getTotalSize( const std::vector<particleGroup> & groups )
{
    int nMax=0;
    for (const auto & group : groups)
    {
        nMax=std::max(group.iEndExtended,nMax);
    }
    
    return nMax+1;
}

int getNParticles( const std::vector<particleGroup> & groups )
{
    int N=0;
    for (const auto & group : groups)
    {
        N+=group.iEnd - group.iStart + 1;
    }

    return N;
}

chain::chain() :
prev(-1) , next(-1) , head(0) , tail(-1)
    {

    }

void pimcConfigurations::setHeadTail( int iChain, int newHead , int newTail )
{
    int oldHead=_chains[iChain].head;
    int oldTail=_chains[iChain].tail;

    int deltaHead = newHead - oldHead;
    int deltaTail = newTail - oldTail;

    if (deltaHead > 0)
    {
        setHead(iChain,newHead);
        setTail(iChain,newTail);
    }
    else
    {
        setTail(iChain,newTail);
        setHead(iChain,newHead);
    }
    
}

void pimcConfigurations::setHead( int iChain, int newHead )
{
    int oldHead=_chains[iChain].head;

    assert(iChain != -1);


    auto & group = getModifiableGroupByChain(iChain);


    int delta = newHead - oldHead;
    
    if ( delta > 0)
    {
        createBeads({oldHead,newHead-1},iChain);
    }
    else
    {
        deleteBeads({newHead,oldHead},iChain);
    }
    
    auto & currentChain=_chains[iChain];
    int oldNext=currentChain.next;

    currentChain.head=newHead;
    currentChain.next=-1;

    assert(currentChain.tail < currentChain.head);

    if ( oldNext  >= 0 )
    {
        setTail(oldNext,_chains[oldNext].tail);
        group.pushHead( iChain);
    }

};


std::list<int> pimcConfigurations::buildPolimerList(int iChain) const
{
    std::list<int> chains;

    assert(iChain != -1);

    
    int iCurrentChain=iChain;
    do
    {
        chains.emplace_back(iCurrentChain);
        iCurrentChain=getChain(iCurrentChain).next ;
    } while (
        (iCurrentChain != -1) and
        (iCurrentChain != iChain)    
        );
    if (iCurrentChain == -1 )
    {
        iCurrentChain=getChain(iChain).prev;
        while(iCurrentChain!= - 1)
        {
            chains.emplace_front(iCurrentChain);
            iCurrentChain=getChain(iCurrentChain).prev;
        }

    }
    return chains;
}

void pimcConfigurations::setTail( int iChain, int newTail )
{
    int oldTail=_chains[iChain].tail;

    assert(iChain != -1 );


    auto & group = getModifiableGroupByChain(iChain);


    int delta = newTail - oldTail;

    if ( delta < 0)
    {
        createBeads({newTail+1,oldTail},iChain);
    }
    else
    {
        deleteBeads({oldTail+1,newTail},iChain);
    }

    auto & currentChain=_chains[iChain];

    int oldPrev=currentChain.prev;

    currentChain.tail=newTail;
    currentChain.prev=-1;

    assert(currentChain.tail < currentChain.head);

    if ( oldPrev  >= 0 )
    {
        setHead(oldPrev,_chains[oldPrev].head);
        group.pushTail( iChain);
    }

};

pimcConfigurations::pimcConfigurations(
    size_t timeSlices, int dimensions, const std::vector<particleGroup> &  particleGroups_) :
    particleGroups(particleGroups_),
    M(timeSlices),
    N(getTotalSize(particleGroups_)), // number of chains(includign padding chains)
    _data(getTotalSize(particleGroups_),dimensions,2*(timeSlices+1) ),// contains a copy for buffer operations
     _tags(getTotalSize(particleGroups_), timeSlices ),
    _nParticles(getNParticles(particleGroups_)),// number of chains without the padding
    mu(0)
{
    _chains.resize(N);

    for(int i=0;i<N;i++)
    {
        _chains[i].tail=-1;
        _chains[i].head=nBeads();
        _chains[i].next=i;
        _chains[i].prev=i;
    }

    _tags.setConstant(0);
    for (const auto & group : particleGroups)
    {
        for (int t=0;t<nBeads();t++)
        {
            for(int i=group.iStart;i<=group.iEndExtended;i++)
            {
                _tags(i,t)=1;
            }
        }
    }


    ensamble = ensamble_t::canonical;

}; 

void pimcConfigurations::setChemicalPotential(const std::vector<Real> & mu_)
{
    mu=mu_;
    setEnsamble(ensamble_t::grandCanonical);
}

void pimcConfigurations::setDifferenceChemicalPotential(const std::vector<Real> & mu_)
{
    if ( getGroups().size() != mu_.size() + 1 )
    {
        throw std::runtime_error("Number of groups should be equal to number of difference chemical potentials minus one.");
    }

    mu.resize(mu_.size() + 1 , 0 );

    for(int i=0;i<mu_.size();i++)
    {
        mu[i+1]= mu_[i];
    }

    setEnsamble(ensamble_t::semiGrandCanonical);

}



void pimcConfigurations::fillHead(int iChain)
{
    const auto & currentChain = _chains[iChain];

    if ( currentChain.next != -1 )
    {

        for (int d=0;d<getDimensions();d++)
        {
            _data(iChain,d,currentChain.head )=_data(currentChain.next,d, 
            _chains[currentChain.next].tail + 1
            );
        }

    }

}

void pimcConfigurations::fillHeads()
{
   for( const auto & group : particleGroups )
   {
       for(int i=group.iStart;i<=group.iEnd;i++)
       {
            fillHead(i);
       }
       
   }    

}

int pimcConfigurations::nParticles() const
{
    int nTot=0;
    for( int i=0;i < particleGroups.size() ; i++ )
    {
      nTot+=nParticles(i);
    }    

    return nTot;
}


int pimcConfigurations::nParticles( int i) const
{
    assert(i < particleGroups.size());
    const auto & group = getGroups()[i];

    
    int N=(group.iEnd - group.iStart + 1);

    if ( group.isOpen() )
    {
        if ( group.heads[0] == group.tails[0] )
        {
            N-=1;
        } 
        else
        {
            N-=2;
        }
    }

    return N;

}

void pimcConfigurations::deleteBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _tags(iChain,t)=0;
    }
}

void pimcConfigurations::createBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _tags(iChain,t)=1;
    }
}

int pimcConfigurations::pushChain( particleGroup & group)
{
    group.iEnd+=1;

    if (group.iEnd > group.iEndExtended )
    {
        throw invalidState("No availible memory left for this group. ");
    }

    if ( _chains[group.iEnd].hasHead() )
    {
        group.pushHead(group.iEnd);
    }

    if ( _chains[group.iEnd].hasTail() )
    {
        group.pushTail(group.iEnd);
    }
    
    return group.iEnd;
}

int pimcConfigurations::pushChain( int iGroup)
{
    auto & group = particleGroups[iGroup];
    return pushChain(group);
}

void pimcConfigurations::removeChain( int iChain)
{
        assert(iChain != -1 );

    auto & group=getModifiableGroupByChain(iChain);

    assert(iChain<= group.iEnd);
    
    // unregister chain
    if ( !(getChain(iChain).hasHead()) or !( getChain(iChain).hasTail()) )
    {
        throw invalidState("Cannot remove chain with no head or tail");
    }

        swap(iChain,group.iEnd);
    
    

    deleteTailFromList(group.iEnd);
    deleteHeadFromList(group.iEnd);

    group.iEnd-=1;
    
}

void pimcConfigurations::removeChains( int iChain1, int iChain2)
{
    assert(iChain1 != -1 );
    assert(iChain2 != -1 );

    auto & group1=getModifiableGroupByChain(iChain1);
    auto & group2=getModifiableGroupByChain(iChain2);

    removeChain(iChain1);

    if (iChain1 == iChain2) 
    {
        return ;
    }


    if (iChain2 !=( group1.iEnd + 1))
    {
        removeChain(iChain2);
    }
    else
    {
        removeChain(iChain1);
    }

}

void pimcConfigurations::deleteHeadFromList(int iChain)
{
    // delete the head from the list of heads
    assert(iChain != -1 );
    auto & group = getModifiableGroupByChain(iChain);

    group.removeHead(iChain);


   
}

void pimcConfigurations::deleteTailFromList(int iChain)
{
        assert(iChain != -1 );
    auto & group = getModifiableGroupByChain(iChain);
    group.removeTail(iChain);

}


void pimcConfigurations::join( int iChainLeft, int iChainRight)
{
    setHead(iChainLeft,_chains[iChainLeft].head);
    setTail(iChainRight,_chains[iChainRight].tail);
    

    
    _chains[iChainLeft].next=iChainRight;
    _chains[iChainRight].prev=iChainLeft;

    deleteHeadFromList(iChainLeft);
    deleteTailFromList(iChainRight);   

}

void pimcConfigurations::swapTails(int iChain1, int iChain2)
{
    auto headChain=getChain(iChain1);
    auto partnerChain = getChain(iChain2);

    setTail(iChain2,headChain.tail);
    setTail(iChain1,partnerChain.tail);

        if ( partnerChain.prev != -1)
        {
            join(partnerChain.prev,headChain.tail);
        }
        if ( headChain.prev != -1)
        {
            join(headChain.prev,headChain.next);
        }
}


std::ostream & pimcConfigurations::operator>>( std::ostream & f) const
{
        const std::string delim = "\t";

        f << "particle"<<" time" << delim << "x" << delim << "y" << delim << "z" <<delim << "tag" << std::endl;
        f << std::setprecision(7);

        for (int t=0;t<=nBeads() ;t++ )
        {
                for (int i=0;i<nChains();i++)
                {
                    
                
                    f<< i << delim << t << delim;

                    #if DIMENSIONS == 3    
                        f <<  _data(i,0,t) << delim;
                        f <<  _data(i,1,t) << delim;
                        f <<  _data(i,1,t) << delim;   
                    #endif

                    #if DIMENSIONS == 1    
                        f <<  _data(i,0,t) << delim;
                        f <<  0 << delim;
                        f << 0 << delim;   
                    #endif

                    #if DIMENSIONS == 2    
                        f <<  _data(i,0,t) << delim;
                        f <<  _data(i,1,t) << delim;
                        f <<  0 << delim;   
                    #endif



                    f << _tags(i,t) << std::endl;
                }
        }

        return f;
}

void pimcConfigurations::save(const std::string & dirname,const std::string & format ) const
{
    if ( ! fs::exists(dirname) ) 
    { 
        fs::create_directory(dirname); // create src folder
    }

    std::ofstream f;

    if ( format == "csv" )
    {
        const std::string delim = "\t";


        f.open(dirname + "/particles.dat");
        //f.write(reinterpret_cast<char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));
        f << "particle"<<" time" << delim << "x" << delim << "y" << delim << "z" <<delim << "mask" << std::endl;
        f << std::setprecision(7);
        
        for (int t=0;t<=nBeads() ;t++ )
        for ( const auto & group : getGroups() )
        {
                for (int i=group.range()[0];i<=group.range()[1];i++)
                {
                
                    f<< i << delim << t << delim;

                    #if DIMENSIONS == 3    
                        f <<  _data(i,0,t) << delim;
                        f <<  _data(i,1,t) << delim;
                        f <<  _data(i,2,t) << delim;   
                    #endif

                    #if DIMENSIONS == 1    
                        f <<  _data(i,0,t) << delim;
                        f <<  0 << delim;
                        f << 0 << delim;   
                    #endif

                    #if DIMENSIONS == 2    
                        f <<  _data(i,0,t) << delim;
                        f <<  _data(i,1,t) << delim;
                        f <<  0 << delim;   
                    #endif



                    f << _tags(i,t) << std::endl;
                }
        }
        f.close();
        f.open(dirname + "/chains.dat");
        f << "chain" << delim << "prev" << delim << "next"<< delim << "head" << delim << "tail"<< std::endl;

         // network information
        for (const auto & group : particleGroups)
        {
            for (int i=group.iStart;i<=group.iEnd;i++)
            {
                auto chain = getChain(i);
                f << i << delim << chain.prev << delim << chain.next << delim<< chain.head << delim << chain.tail << std::endl;
            }
        }

    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.write(reinterpret_cast<const char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));

        throw missingImplementation("Binary format not yet supported.");
    }
    else if (format == "pdb")
    {
        int k=0;
        f.open(dirname + "/particles.pdb");

        f << std::fixed << std::setprecision(3) << std::right ;
        for (int t=0;t<=nBeads() ;t++ )
        {
             for (const auto & group : getGroups())
                for (int i=group.iStart;i<=group.iEndExtended;i++)
                {
                    
                    bool isVisible = (i <= group.iEnd) ? 1 : 0;
                    
                    isVisible = isVisible and (bool)(std::max(_tags(i,t) , _tags(i,t-1) ));

                    Real bogus=1e+3;

                    Real x = 0;
                    Real y = 0;
                    Real z = 0;
                    
                    if (isVisible)
                    {
                        # if DIMENSIONS == 3
                        x = _data(i,0,t);
                        y = _data(i,1,t);
                        z = _data(i,2,t);
                        #endif

                        # if DIMENSIONS == 2
                        x = _data(i,0,t);
                        y = _data(i,1,t);
                        #endif

                        # if DIMENSIONS == 1
                        x = _data(i,0,t);
                        #endif
                    }
                    else
                    {
                        # if DIMENSIONS == 3
                        x = bogus;
                        y = bogus;
                        z = bogus;
                        #endif

                        # if DIMENSIONS == 2
                        x = bogus;
                        y = bogus;
                        #endif

                        # if DIMENSIONS == 1
                        x = bogus;
                        #endif
                    }
                    


                    f<< "HETATM"
                    <<  std::setw(5) << nChains()*t + i << " C    VAL A" << std::setw(4) <<  t << "    "
                    
                    << std::setw(8) << x
                    << std::setw(8) << y
                    << std::setw(8) << z
                    
                    
                    
                    << std::setw(6) << std::setprecision(2) << 1.0
                    << std::setw(6) << std::setprecision(2) << i 
                    << "     C" << std::endl;

                    k++;    
                }
                
        }

        for (int i=0;i<nChains();i++)
            for (int t=0;t<nBeads() ;t++ )
            {
                {
                    f<< "CONECT"
                    <<  std::setw(5) << nChains()*t + i 
                    <<  std::setw(5) << nChains()*(t+1) + i 
                    << std::endl;
                    k++;    
                }
                
            }



        f.close();
    }


    nlohmann::json j;
    j["nChains"]= nChains();
    j["dimensions"]=getDimensions();
    j["timeSlices"]=nBeads();
    j["format"]=format;

    std::vector<nlohmann::json> jGroups;
    for( const auto & group  : particleGroups)
    {
        nlohmann::json jGroup;
        jGroup["iStart"]=group.iStart;
        jGroup["iEnd"]=group.iEnd;
        jGroup["iEndExtended"]=group.iEndExtended;
        jGroups.push_back(jGroup);
    }
    
    j["groups"]=jGroups;
    f.open(dirname + "/description.json");
    f << j ;
    f.close();

}

void pimcConfigurations::load(const std::string & dirname)
{

    if ( ! fs::exists(dirname) ) 
    { 
        throw invalidInput("Folder " + dirname + " does not exist.");
    }

    std::ifstream f;
    f.open(dirname + "/description.dat");
    nlohmann::json j;
    
    f >> j;

    f.close();


    N = j["nChains"].get<int>();
    M=j["timeSlices"].get<int>();

    _data.resize(N,getDimensions(),M);
    _tags.resize(nChains(),nBeads()) ;
    
    _tags.setConstant(0);


    auto format = j["format"].get<std::string>();

    if (format == "csv")
    {
        std::string dummy;
        const int nFields = 5;
        f.open(dirname + "/particles.dat");
        for(int i=0;i<nFields ; i++)
        {
            f >> dummy;
        }

        //f.read(reinterpret_cast<char*>(_data.data()),N*getDimensions()*M*sizeof(double));
        int ii=0, tt=0;
        _data.setConstant(0);

        while(!f.eof() )
                {
                    f >> ii;
                    f >> tt;
                    for (int d=0;d<getDimensions();d++)
                    {
                        f >> _data(ii,d,tt) ;
                        _tags(ii,tt)=1;
                    }
                }
    

    f.close();
    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.read(reinterpret_cast<char*>(_data.data()),N*getDimensions()*M*sizeof(double));
    }

    throw missingImplementation("Chain info not saved");
    
    // set up the groups
    particleGroups={};
    
    for (const auto & group : j["groups"])
    {
       particleGroups.push_back(
           particleGroup(
           {group["iStart"].get<int>(),group["iEnd"].get<int>(),group["iEndExtended"].get<int>(),1.0
           }
           ));
    }

    

}

void pimcConfigurations::copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        Eigen::Tensor<Real,3> & dataTo,
       int timeOffsetTo, int particleOffsetTo )
        {
            const auto & dataFrom = confFrom.dataTensor();

            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        dataTo(ii,d,tt)=dataFrom(i,d,t);
                    }
        }

        }


void pimcConfigurations::translateData(const std::array<int,2> & timeRange, const std::array<int,2> & particleRange,const std::array<Real,getDimensions()> & delta)
{

    for(int t=timeRange[0];t<=timeRange[1];t++)
    {
        
            for(int i=particleRange[0];i<=particleRange[1];i++)
            {
                for(int d=0;d<getDimensions();d++)
                {
                    _data(i,d,t)+=delta[d];
                }
                
            }
    }

}


void pimcConfigurations::copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset) const
        {

            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        buffer(tt,d)=_data(iParticle,d,t);
                    }
            }
        
        }

void pimcConfigurations::copyDataFromBuffer( const Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset)
        {

            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        _data(iParticle,d,t)=buffer(tt,d);
                    }
            }
        }

void pimcConfigurations::swapData( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo)
        {
             auto & dataFrom = confFrom.dataTensor();
            auto & dataTo = confTo.dataTensor();
            

            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        std::swap(dataTo(ii,d,tt),dataFrom(i,d,t));
                    }
            }
        }


void pimcConfigurations::swapTags( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo)
        {
            auto & dataFrom = confFrom._tags;
            auto & dataTo = confTo._tags;

            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        std::swap(dataTo(ii,tt),dataFrom(i,t));
                    }
            }
        }



void pimcConfigurations::swap(int particleA, int particleB)
{

    if (particleA == particleB)
    {

        return;
    }

    // relink chains under the swap
    auto chainA=_chains[particleA];
    auto chainB=_chains[particleB];

    assert( particleA != -1 );
    assert( particleB != -1 );


    auto & groupA = getModifiableGroupByChain(particleA);
    auto & groupB = getModifiableGroupByChain(particleB);

    assert(groupA == groupB); // only swap between particles of the same set are allowed

    if (chainA.prev != -1)
    {
        _chains[chainA.prev].next=particleB; 
    }
    if (chainA.next != -1)
    {
        _chains[chainA.next].prev=particleB; 
    }

    if (chainB.prev != -1)
    {
        _chains[chainB.prev].next=particleA; 
    }
    if (chainB.next != -1)
    {
        _chains[chainB.next].prev=particleA; 
    }


    // swap the chain info
    std::swap(_chains[particleA],_chains[particleB]);

    auto & _heads = groupA.heads;
    auto & _tails = groupA.tails;


    // swap A and B indices in the head list
    std::replace(_heads.begin(),_heads.end(),particleA,-1);
    std::replace(_heads.begin(),_heads.end(),particleB,particleA);
    std::replace(_heads.begin(),_heads.end(),-1,particleB);


    // swap A and B indices in the tail list
    std::replace(_tails.begin(),_tails.end(),particleA,-1);
    std::replace(_tails.begin(),_tails.end(),particleB,particleA);
    std::replace(_tails.begin(),_tails.end(),-1,particleB);

    // swap data

    swapData(particleA,particleB);

    if ( ensamble != ensamble_t::canonical)
    {
        swapTags(particleA,particleB);
    }


}



int configurationsSampler::sampleChain(const configurations_t & confs,randomGenerator_t & randG)
{
    // sample a chain with probability 1/ N_particles 
    int iParticle = uniformRealNumber(randG)*confs.nParticles();
    int k=0;
    int iChain=-1;

    for(const auto & group : confs.getGroups() )
    {
        k+=group.size();
        
        if (k> iParticle)
        {
            iChain = group.iEnd + 1 - (k-iParticle);
        }
    }

    return iChain;
}


int configurationsSampler::sampleChain(const configurations_t & confs, int iGroup,randomGenerator_t & randG)
{
    // sample a chain from a group with probability 1/ (# particles in the group) 

    const auto & group = confs.getGroups()[iGroup];

    if (group.size() == 0)
    {
        return -1;
    }

    
    int iParticle = std::floor(uniformRealNumber(randG)*group.size());
    int k=0;
    int iChain=-1;

    iChain = group.iStart + iParticle;

    return iChain;
}


void configurationsSampler::sampleFreeParticlePosition(
    std::array<Real,getDimensions()> & x,const std::array<Real,getDimensions()> & mean,Real tau,randomGenerator_t & randG,Real mass
){
    Real var = 2 * D * tau / mass;
    for(int d=0;d<getDimensions();d++)
    {
        x[d]=mean[d] + normal(randG)*std::sqrt(var);       
    }
}

int configurationsSampler::sampleGroup(const configurations_t & confs,randomGenerator_t & randG)
{
     return  std::floor(uniformRealNumber(randG)*confs.getGroups().size());
}


void pimcConfigurations::saveHDF5(const std::string & filename)
{
    hid_t status;
    hid_t file = H5Fcreate (filename.c_str() , H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    hsize_t dims[3];

    for(int d=0;d<3;d++)
    {
        dims[2 - d]=dataTensor().dimensions()[d];
    }

    hsize_t scalarDim {1} ;
     hid_t scalarSpace = H5Screate_simple (1 , &scalarDim, NULL);

    int M = nBeads();
    auto nBeads_id= H5Acreate (file,"nBeads", H5T_NATIVE_INT, scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite ( nBeads_id , H5T_NATIVE_INT, & M );
    H5Aclose(nBeads_id);
    


    int ensamble = getEnsamble();

    auto ensamble_id= H5Acreate (file,"ensamble", H5T_NATIVE_INT, scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite ( ensamble_id , H5T_NATIVE_INT, & ensamble );
    H5Aclose( ensamble_id );

    setEnsamble( (ensamble_t)ensamble);
    if (getEnsamble() == ensamble_t::grandCanonical)
    {
        hsize_t muDims { getGroups().size() } ;
        hid_t muSpace = H5Screate_simple ( 1 , &muDims, NULL);
        auto mu_id= H5Acreate (file,"mu", H5T_NATIVE_DOUBLE, muSpace, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Awrite ( mu_id , H5T_NATIVE_DOUBLE, mu.data() );
        H5Aclose( mu_id );
        H5Sclose(muSpace);

    }


    
    hid_t tensorSpace  = H5Screate_simple (3, dims, NULL);

    
    for (int iGroup=0;iGroup<getGroups().size() ; iGroup++  )
    {
        const auto & group = getGroups()[iGroup];

        std::string groupName = "/set" + std::to_string(iGroup);
        hid_t group_id = H5Gcreate(file, groupName.c_str() , H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        auto nGroup = group.iEnd - group.iStart + 1;

        hsize_t dimsGroup[3];
        dimsGroup[2]=nGroup;
        dimsGroup[1]=DIMENSIONS;
        dimsGroup[0]=nBeads();

        hsize_t offset[3], count[3];

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = group.iStart;

        count[0]=nBeads();
        count[1]=DIMENSIONS;
        count[2] = nGroup;

        status = H5Sselect_hyperslab( tensorSpace, H5S_SELECT_SET, offset, NULL, count, NULL);

        hid_t groupSpace  = H5Screate_simple (3, dimsGroup, NULL);

        hid_t dset = H5Dcreate(group_id, "particles", H5T_NATIVE_DOUBLE, groupSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

        status = H5Dwrite (dset, H5T_NATIVE_DOUBLE, tensorSpace, groupSpace,H5P_DEFAULT,  data() );

        hid_t iStart_id= H5Acreate (group_id,"iStart", H5T_NATIVE_INT, scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
        hid_t iEnd_id= H5Acreate (group_id,"iEnd", H5T_NATIVE_INT, scalarSpace, H5P_DEFAULT, H5P_DEFAULT);
        hid_t iEndExtended_id= H5Acreate (group_id,"iEndExtended", H5T_NATIVE_INT, scalarSpace, H5P_DEFAULT, H5P_DEFAULT);

        status = H5Awrite (iStart_id, H5T_NATIVE_INT, &group.iStart );
        status = H5Awrite (iEnd_id, H5T_NATIVE_INT, &group.iEnd );
        status = H5Awrite (iEndExtended_id, H5T_NATIVE_INT, &group.iEndExtended );

        std::vector<int> tails(nGroup,0);
        std::vector<int> heads(nGroup,0);
        std::vector<int> nexts(nGroup,0);
        std::vector<int> prevs(nGroup,0);


        for(int i=0;i<nGroup;i++)
        {
            tails[i]=getChain(i+group.iStart).tail;
            heads[i]=getChain(i+group.iStart).head;
            nexts[i]=getChain(i + group.iStart).next;
            prevs[i]=getChain(i+group.iStart).prev;
        }

        hsize_t particleDims{nGroup};
        auto particleSpace  = H5Screate_simple(1,&particleDims, NULL);


        hid_t heads_id = H5Dcreate(group_id, "heads", H5T_NATIVE_INT, particleSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (heads_id, H5T_NATIVE_INT, particleSpace, particleSpace,H5P_DEFAULT,  heads.data() );
        H5Dclose(heads_id);

        hid_t tails_id = H5Dcreate(group_id, "tails", H5T_NATIVE_INT, particleSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite (tails_id, H5T_NATIVE_INT, particleSpace, particleSpace,H5P_DEFAULT,  tails.data() );
        H5Dclose(tails_id);

        hid_t nexts_id = H5Dcreate(group_id, "nexts", H5T_NATIVE_INT, particleSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite ( nexts_id , H5T_NATIVE_INT, particleSpace, particleSpace,H5P_DEFAULT,  nexts.data() );
        H5Dclose(nexts_id);

        hid_t prevs_id = H5Dcreate(group_id, "prevs", H5T_NATIVE_INT, particleSpace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite ( prevs_id , H5T_NATIVE_INT, particleSpace, particleSpace,H5P_DEFAULT,  prevs.data() );
        H5Dclose(prevs_id);

        status = H5Sclose (particleSpace);
        
        status = H5Aclose (iStart_id);
        status = H5Aclose (iEnd_id);
        status = H5Aclose (iEndExtended_id);


        status = H5Dclose (dset);
        status = H5Sclose (groupSpace);


        status = H5Gclose ( group_id);

    }



    status = H5Sclose (tensorSpace);
    status = H5Sclose ( scalarSpace);
        
    status = H5Fclose (file);
}


herr_t listGroups(hid_t loc_id, const char *name, const H5L_info_t *info, void * data)
    {
        auto vec= (std::vector<std::string > *) (data);

        vec->push_back( std::string( name ));

        return 0;

        
    }
    


pimcConfigurations pimcConfigurations::loadHDF5(const std::string & filename)
{
    hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

    
    std::vector<std::string> groupNames;
    

    std::vector<particleGroup> groups;


    
    auto status=H5Literate (file_id, H5_INDEX_NAME, H5_ITER_INC, NULL, listGroups ,(void *) &groupNames);


    int nBeads;
    auto nBeads_id =H5Aopen( file_id, "nBeads", H5P_DEFAULT );
    H5Aread(nBeads_id, H5T_NATIVE_INT, &nBeads);
    H5Aclose(nBeads_id);


    for ( const auto &  groupName : groupNames)
    {
        auto group_id=H5Gopen( file_id, groupName.c_str() ,H5P_DEFAULT);

        auto iStart_id =H5Aopen( group_id, "iStart", H5P_DEFAULT );
        auto iEnd_id =H5Aopen( group_id, "iEnd", H5P_DEFAULT );
        auto iEndExtended_id =H5Aopen( group_id, "iEndExtended", H5P_DEFAULT );

        int iStart,iEnd,iEndExtended;
        auto status = H5Aread(iStart_id, H5T_NATIVE_INT, &iStart);
        status = H5Aread(iEnd_id, H5T_NATIVE_INT, &iEnd);
        status = H5Aread(iEndExtended_id, H5T_NATIVE_INT, &iEndExtended);

        groups.emplace_back( iStart,iEnd,iEndExtended   );

        H5Aclose(iStart_id);
        H5Aclose(iEnd_id);
        H5Aclose(iEndExtended_id);

        /* auto dataSet_id=H5Dopen(group_id,"particles", H5P_DEFAULT );
        auto space_id=H5Dget_space( dataSet_id );

        hsize_t dims[3];

        int rank= H5Sget_simple_extent_dims( space_id, dims, NULL ); 

        assert(rank==3);
        

        H5Dclose( dataSet_id);
        H5Dclose(space_id);
        */

       H5Gclose(group_id);
    }

    pimc::pimcConfigurations configurations(nBeads,DIMENSIONS,groups);

     int ensamble;
    auto ensamble_id =H5Aopen( file_id, "ensamble", H5P_DEFAULT );
    H5Aread(ensamble_id, H5T_NATIVE_INT, &ensamble);
    H5Aclose(ensamble_id);
    
     if (ensamble == ensamble_t::grandCanonical)
    {
        std::vector<Real> mu(groups.size(),0);
        auto mu_id =H5Aopen( file_id, "mu", H5P_DEFAULT );
        H5Aread(mu_id, H5T_NATIVE_DOUBLE, mu.data() );
        H5Aclose(mu_id);

        configurations.setChemicalPotential(mu);
    }


    

    hsize_t dims[3];
    for(int d=0;d<3;d++)
    {
        dims[2 - d]=configurations.dataTensor().dimensions()[d];
    }
    hid_t tensorSpace  = H5Screate_simple (3, dims, NULL);



    for(int iGroup=0;iGroup<groups.size();iGroup++)
    {
        const auto & group = configurations.getGroups()[iGroup];

        auto group_id=H5Gopen( file_id,groupNames[iGroup].c_str() ,H5P_DEFAULT);

        auto dataSet_id=H5Dopen(group_id,"particles", H5P_DEFAULT );
        auto space_id=H5Dget_space( dataSet_id );


      
        hsize_t offset[3], count[3];
        auto nGroup = group.iEnd - group.iStart + 1;

        offset[0] = 0;
        offset[1] = 0;
        offset[2] = group.iStart;
        count[0]=configurations.nBeads();
        count[1]=DIMENSIONS;
        count[2] = nGroup;
        status = H5Sselect_hyperslab( tensorSpace, H5S_SELECT_SET, offset, NULL, count, NULL );

        auto status=H5Dread( dataSet_id, H5T_NATIVE_DOUBLE, tensorSpace, space_id, H5P_DEFAULT, configurations.data() );


    /*     hsize_t dims[3];
        int rank= H5Sget_simple_extent_dims( space_id, dims, NULL ); 
        assert(rank==3); */

       
        std::vector<int> heads(nGroup,0);
        std::vector<int> tails(nGroup,0);
        std::vector<int> nexts(nGroup,0);
        std::vector<int> prevs(nGroup,0);

        auto heads_id=H5Dopen(group_id,"heads", H5P_DEFAULT );
        auto particleSpace = H5Dget_space( heads_id );
        status=H5Dread( heads_id, H5T_NATIVE_INT, particleSpace,particleSpace, H5P_DEFAULT, heads.data() );
        H5Dclose(heads_id);

        auto tails_id=H5Dopen(group_id,"tails", H5P_DEFAULT );
        status=H5Dread( tails_id, H5T_NATIVE_INT, particleSpace,particleSpace, H5P_DEFAULT, tails.data() );
        H5Dclose(tails_id);

        auto nexts_id=H5Dopen(group_id,"nexts", H5P_DEFAULT );
        status=H5Dread( nexts_id, H5T_NATIVE_INT, particleSpace,particleSpace, H5P_DEFAULT, nexts.data() );
        H5Dclose(nexts_id);

        auto prevs_id=H5Dopen(group_id,"prevs", H5P_DEFAULT );
        status=H5Dread( prevs_id, H5T_NATIVE_INT, particleSpace,particleSpace, H5P_DEFAULT, prevs.data() );
        H5Dclose(prevs_id);

        for(int i=0;i<nGroup;i++)
        {
            int iChain=i + group.iStart;
            configurations.setHeadTail(iChain,heads[i],tails[i]);

        }

        for(int i=0;i<nGroup;i++)
        {
            int iChain=i + group.iStart;

            if (nexts[i] >=0 )
            {
                configurations.join(iChain,nexts[i]);
                
            }
        }
        





        H5Sclose(particleSpace);
        H5Dclose( dataSet_id);
        
        H5Sclose(space_id);
        H5Gclose(group_id);


    }

    H5Sclose(tensorSpace);
    H5Fclose(file_id);

    configurations.fillHeads();


    return configurations;

}

int nParticlesOnClose(const pimcConfigurations & configurations, int set)
{
    int nParticles=configurations.nParticles(set);
    int M = configurations.nBeads();
    if ( configurations.getGroups()[set].isOpen() )
    {
        const auto & group = configurations.getGroups()[set];
        int iChainHead=group.heads[0];
        int iChainTail=group.tails[0];
        int tHead = configurations.getChain(iChainHead).head;
        int tTail = configurations.getChain(iChainTail).tail;
        
        if (iChainTail == iChainHead)
        {
            nParticles+=1;
            if (tHead - (tTail + 1 ) >=M)
            {
                nParticles+=1;
            }
        }
        else
        {
            if (tHead >= (tTail + 1) )
            {
                nParticles+=2;
            }
            else
            {
                nParticles+=1;
            }
        }

    }

    
    return nParticles;
    
}



void pimcConfigurations::setRandom( const std::array<Real,DIMENSIONS> & lBox,randomGenerator_t & randG)
{
    std::uniform_real_distribution<double> uniformDistribution(0.0,1.0);

    auto & data=dataTensor();

    for (int t=0;t<data.dimensions()[2];t++)
        for (int i=0;i<data.dimensions()[0];i++)
            for  (int d=0;d<DIMENSIONS;d++)
            {
                data(i,d,t)=(uniformDistribution(randG)-0.5 )*lBox[d];
            }
    
    fillHeads();

}
    
};



