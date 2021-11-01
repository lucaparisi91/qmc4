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
    int rank = 1;
    const auto & data = dataTensor();

    size_t dims[rank];

    dims[0]=data.dimensions()[0]*data.dimensions()[1] * data.dimensions()[2] ;


    const auto & groups =getGroups();
    int chunkSize=6;
    std::vector<int> groupData(groups.size() * chunkSize );

    int groupDims[1]={groups.size()*chunkSize};

    std::vector<double> masses;
    masses.resize(groups.size());

    for(int i=0;i<groups.size();i++)
    {
        groupData[i*chunkSize]=groups[i].iStart;
        groupData[i*chunkSize+1]=groups[i].iEnd;
        groupData[i*chunkSize+2]=groups[i].iEndExtended;
        
        if ( groups[i].isOpen() )
        {
            groupData[i*chunkSize+4]=groups[i].heads[0];
            groupData[i*chunkSize+5]=groups[i].tails[0];
        }
        else
        {
            groupData[i*chunkSize+4]=-1;
            groupData[i*chunkSize+5]=-1;
        }
        
        
        masses[i]=groups[i].mass;
    }


    std::vector<int> chainData;
    int chainChunkSize=4;

    chainData.resize(_chains.size()*chainChunkSize);
    for(int i=0;i<_chains.size();i++)
    {
        chainData[chainChunkSize*i ]=_chains[i].prev;
        chainData[chainChunkSize*i + 1 ]=_chains[i].next;
        chainData[chainChunkSize*i + 2]=_chains[i].head;
        chainData[chainChunkSize*i + 3]=_chains[i].tail;
    }

    hdf5IO ioInterface(filename, std::ios::out  );

    ioInterface.write(data.data(),"configurations",& dims[0],rank);
    ioInterface.write(groupData,"groupings");
    ioInterface.write(chainData,"chains");

    ioInterface.write(masses,"mass");
    ioInterface.annotate("nBeads",M,"configurations");
    ioInterface.close();

}

pimcConfigurations pimcConfigurations::loadHDF5(const std::string & filename)
{
    hdf5IO ioInterface2(filename,std::ios::in | std::ios::out );

    auto groupData2=ioInterface2.get<std::vector<int> >("groupings");
    auto masses2=ioInterface2.get<std::vector<double> >("mass");
    int M2[1];
    ioInterface2.readNote("nBeads",M2,"configurations");

    std::cout << "Open file " << filename << std::endl;

    std::vector<pimc::particleGroup> groups2;
    int chunkSize=6;
    for(int i=0;i<masses2.size();i++)
    {
        pimc::particleGroup currentGroup ( 
            groupData2[i*chunkSize], 
            groupData2[i*chunkSize + 1],
            groupData2[i*chunkSize + 2],
            masses2[i]
        );  

        int iHead =groupData2[i*chunkSize + 4] ;
        int iTail = groupData2[i*chunkSize + 5] ;

        if ( ( iTail != -1 ) or ( iHead != -1 ) )
        {
            currentGroup.pushHead(iHead);
            currentGroup.pushTail(iTail);
        } 
        else 
        {
            // check that both tail and head are missing
            assert(iTail == -1);
            assert(iHead == -1);
        }

        groups2.push_back(currentGroup);
    }

    
    pimc::pimcConfigurations configurations2(M2[0] , getDimensions() , groups2); 
    //pimc::pimcConfigurations configurations2(50 , getDimensions() , {{0,999,999,1}});

    ioInterface2.read(configurations2.data(), "configurations"  );


    auto chainData=ioInterface2.get<std::vector<int> >("chains");

    int chainChunkSize=4;
    for(int i=0;i<configurations2._chains.size() ;i++ )
    {
        configurations2._chains[i].prev=chainData[i*chainChunkSize];
        configurations2._chains[i].next=chainData[i*chainChunkSize+1];
        configurations2._chains[i].head=chainData[i*chainChunkSize+2];
        configurations2._chains[i].tail=chainData[i*chainChunkSize+3];
    }

    ioInterface2.close(); 
    return configurations2;
}

int nParticlesOnClose(const pimcConfigurations & configurations, int set)
{
    int nParticles=configurations.nParticles();
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



};



