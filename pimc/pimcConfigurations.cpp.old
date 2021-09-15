#include "pimcConfigurations.h"
#include "qmcExceptions.h"
#include <filesystem>
#include <sstream>
#include <algorithm>

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

pimcConfigurations::pimcConfigurations(
    size_t timeSlices, int dimensions, const std::vector<particleGroup> &  particleGroups_) :
    particleGroups(particleGroups_),
    M(timeSlices),N(getTotalSize(particleGroups_)),
    _data(N,dimensions,timeSlices),
     _mask( timeSlices,N),
    _nParticles(getNParticles(particleGroups_)),
    emptyGroup{0,-1,-1,0.}
{
    for (const auto & group : particleGroups)
    {
        for(int t=0;t<nBeads();t++)
        {
             for(int i=group.iEnd+1;i<=group.iEndExtended;i++)
             {
                 _mask(t,i)=0;
             }
            
        }
    }

}; 

void pimcConfigurations::deleteBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=0;
    }
}

void pimcConfigurations::createBeads(   std::array<int,2> timeRange, int iChain )
{
    for (int t=timeRange[0] ; t<= timeRange[1] ;t++   )
    {
        _mask(t,iChain)=1;
    }
}

int pimcConfigurations::open(   int time, int iChain ,bool copyBack)
{
    if ( _worms.size() > 0 )
    {
        throw missingImplementation("More then on worm is not supported.");

    }
    auto & group = getGroupByChain(iChain);
    
    if(copyBack)
    {
        for(int t=0;t<=time;t++)
        {
            for(int d=0;d<getDimensions();d++)
            {
                _data(group.iEnd+1,d,t)=_data(iChain,d,t);
            }
        }
    }
    
    group.iEnd+=1;
    int iChain2 = group.iEnd; 

    deleteBeads({0,time-1},iChain);
    createBeads({time+1,nBeads()-1},iChain);

    createBeads({0,time},iChain2);
    deleteBeads({time+1,nBeads()-1},iChain2);

    _worms.push_back({iChain2,iChain,time,time});
    return _worms.size()-1;
}

void pimcConfigurations::close( int iWorm, bool copyBack )
{
    const auto & worm = _worms[iWorm];

    auto & group = getGroupByChain(worm.iChainHead);

    // copy the data from the head to the tail
    if (copyBack)
    {
        copyData({0,worm.iHead}, worm.iChainHead,worm.iChainTail);
    }


    createBeads({0,nBeads()-1},worm.iChainTail );
    
    // remove the worm from the proper group

    if (worm.iChainHead != group.iEnd)
    {
        swap(worm.iChainHead,group.iEnd);
    }


    group.iEnd-=1;

    if (group.iEnd<group.iStart)
    {
        throw invalidInput("iEnd is lower than iStart. ");
    };
    
    if (_worms.size() > 1 )
    {
        throw missingImplementation("Multiple worms are not supported");
    };

    _worms.resize(0);
}


void pimcConfigurations::save(const std::string & dirname,const std::string & format )
{

    if ( ! fs::exists(dirname) ) 
    { 
        fs::create_directory(dirname); // create src folder
    }

    std::ofstream f;

    if ( format == "csv" )
    {
    
        f.open(dirname + "/particles.dat");
        //f.write(reinterpret_cast<char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));
        f << "particle time x y z" << std::endl;
        f << std::setprecision(7);

        for (int t=0;t<nBeads() ;t++ )
        {
                for (int i=0;i<nChains();i++)
                {
                    if (_mask(t,i) == 1)
                        {
                            f<< i << " " << t << " ";
                            for (int d=0;d<getDimensions();d++)
                            {
                                
                                f <<  _data(i,d,t) << " ";
                            }
                            f<< std::endl;
                        }
                    
                }
        }
        f.close();
    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.write(reinterpret_cast<char*>(_data.data()), nBeads()*nChains()*getDimensions()*sizeof(double));
    }


    nlohmann::json j;

    j["nChains"]= nChains();
    j["dimensions"]=getDimensions();
    j["timeSlices"]=nBeads();
    j["format"]=format;
    

    std::vector<nlohmann::json> wormObjects;

    for(int i=0;i<_worms.size();i++)
    {
        nlohmann::json jWorm;
        jWorm["iChainHead"]=_worms[i].iChainHead;
        jWorm["iChainTail"]=_worms[i].iChainTail;
        jWorm["iHead"]=_worms[i].iHead;
        jWorm["iTail"]=_worms[i].iTail;
        wormObjects.push_back(jWorm);
    }

    j["worms"]=wormObjects;

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
    f.open(dirname + "/description.dat",std::ios::binary);
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
    _mask=mask(nBeads(),nChains()) ;
    _mask.setConstant(0);

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
                        _mask(tt,ii)=1;
                    }
                }
    

    f.close();
    }
    else if (format == "binary")
    {
        f.open(dirname + "/particles.dat",std::ios::binary);
        f.read(reinterpret_cast<char*>(_data.data()),N*getDimensions()*M*sizeof(double));
    }

    for (const auto & wormJ : j["worms"])
    {
       _worms.push_back( {
           wormJ["iChainHead"].get<int>() ,
           wormJ["iChainTail"].get<int>(),
           wormJ["iHead"].get<int>(),
           wormJ["iTail"].get<int>()
       });
    }


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


void pimcConfigurations::copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset) const
        {
            
            for(int t=timeRange[0], tt=timeOffset;t<=timeRange[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                    {
                        buffer(tt,d)=_data(iParticle,d,t);
                    }
        }
        
        }


void pimcConfigurations::copyDataFromBuffer( const Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset)
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
            auto & maskFrom = confFrom.getMask();
            auto & maskTo = confTo.getMask();


            for(int t=timeRangeFrom[0], tt=timeOffsetTo;t<=timeRangeFrom[1];t++ & tt++)
            {
                for(int d=0;d<getDimensions();d++)
                   for(int i=particleRangeFrom[0], ii=particleOffsetTo;i<=particleRangeFrom[1];i++ & ii++)
                    {
                        std::swap(dataTo(ii,d,tt),dataFrom(i,d,t));
                    }
            }

        }



void pimcConfigurations::updateMask(const worm & wormOld,const worm & wormNew)
{

    // unmask the chains previously covered by worms
    createBeads({0,nBeads()-1},wormOld.iChainHead);
    createBeads({0,nBeads()-1},wormOld.iChainTail);

    // create the masking for the new worm
    createBeads({0,wormNew.iHead},wormNew.iChainHead);
    createBeads({wormNew.iTail,nBeads()-1},wormNew.iChainTail);
    deleteBeads({0,wormNew.iTail-1},wormNew.iChainTail);
    deleteBeads({wormNew.iChainHead+1,nBeads()-1},wormNew.iChainHead);

}

void pimcConfigurations::swap(int particleA, int particleB)
{
    swapData(particleA,particleB);

    // if a worm containts the chain A or B swap the chain
    // Update the mask after updating the worm configuration
    for (int i=0;i<_worms.size();i++)
    {
        auto worm= _worms[i];

        if ( particleA == worm.iChainHead)
        {
            worm.iChainHead=particleB;
        }
        else if (particleA == worm.iChainTail)
        {
            worm.iChainTail=particleB;
        }

        if ( particleB == worm.iChainHead)
        {
            worm.iChainHead=particleA;
        }
        else if (particleB == worm.iChainTail)
        {
            worm.iChainTail=particleA;
        }

        updateWorm(i,worm);

    }

}



};

