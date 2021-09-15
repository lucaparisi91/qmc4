
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS

#include <vector>
#include "../src/traits.h"
#include "../src/tools.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "qmcExceptions.h"

namespace pimc
{
    
class mask
{
public:
    mask(){}
    mask( int nBeads, int nChains) : _mask(nChains,nBeads) {_mask.setConstant(1);}

    void setConstant(int value){_mask.setConstant(value);}


    auto & operator()(int i,int j) {return _mask(j,i);}
    const auto &  operator()(int i,int j) const {return _mask(j,i);}
private:
    Eigen::Tensor<int, 2> _mask;
};

struct particleGroup
{
    
    bool contains(int iParticle) const {return (iParticle>= iStart) and (iParticle<=iEnd);}
    int iStart; // start of the particle group
    int iEnd; // end of the active group
    int iEndExtended; // extended memory for additional particles
    Real mass;

    auto size() const {return iEnd - iStart + 1;}
    
};

struct worm
{
    int iChainHead;
    int iChainTail;
    int iHead;
    int iTail;
};

    class pimcConfigurations
    {
        public:
        using configurationsStorage_t =  Eigen::Tensor<Real, 3> ;

        pimcConfigurations() : pimcConfigurations(0,getDimensions(), {} 
        ) {}

        pimcConfigurations(size_t timeSlices, int dimensions, const std::vector<particleGroup> & particleGroups_);
        bool isOpen() {return _worms.size()> 0 ;}

        auto &  dataTensor() {return _data;}
         const auto &  dataTensor() const {return _data;}


        auto data() {return _data.data();}
        
        auto nChains() const {return N;} 
        
        auto nBeads() const {return M; } 

        const auto & getMask() const {return _mask;}

        const auto & getGroup(int iChain) const
        {
            for (int i=0;i<particleGroups.size();i++)
            {
                if ( particleGroups[i].contains(iChain) )
                {
                    return particleGroups[i];
                }
            }

            throw invalidInput("missing group for particle " + std::to_string(iChain));
        }
        
        int open( int time, int iChain,bool copyBack=true); // adds a worm

        void close(int iWorm,bool copyBack=true); // closes the ith worm (not the the ith chain)

        auto & worms() {return _worms;}

        void updateWorm(int iWorm,const worm &newWorm)
        {
           
            updateMask(_worms[iWorm],newWorm);
             _worms[iWorm]=newWorm;
        };


        void save(const std::string & directoryName,const std::string & format="csv");

        void load(const std::string & directoryName);

        int nParticles() {return _nParticles;}

        const auto & getGroups() const {return particleGroups;}


        static void copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom, const std::array<int,2> & particleRangeFrom ,
          Eigen::Tensor<Real,3> & dataTo, int timeOffsetTo, int particleOffestTo );
  


          static void copyData(const pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom, const std::array<int,2> & particleRangeFrom ,
         pimcConfigurations & confTo, int timeOffsetTo, int     particleOffestTo )
         {
             copyData(confFrom,timeRangeFrom,particleRangeFrom,confTo.dataTensor(),timeOffsetTo,particleOffestTo);

         }


         void copyData(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::copyData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }


         void copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset=0) const;

         void copyDataToBuffer( Eigen::Tensor<Real,2> & buffer, std::array<std::array<int,2>,2> & timeRanges, int iParticle ,int timeOffset=0) 
         {
             copyDataToBuffer(buffer,timeRanges[0],iParticle,timeOffset);
             timeOffset+=timeRanges[0][1]-timeRanges[0][0] + 1;
             copyDataToBuffer(buffer,timeRanges[1],iParticle,timeOffset);
         }

         void copyDataFromBuffer(const  Eigen::Tensor<Real,2> & buffer, std::array<int,2> & timeRange, int iParticle ,int timeOffset=0);

        void swapData( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo );

       void swapData(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::swapData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }

        void swapData(int iParticleFrom, int iParticleTo)
         {
             swapData({0,nBeads()},iParticleFrom,iParticleTo);
         }


        void swap(int iParticleFrom, int iParticleTo);


        particleGroup & getGroupByChain(int iChain) {
            for (auto & group : particleGroups )
            {
                if (group.contains(iChain) )
                {
                    return  group;       
                }
            }
            throw invalidInput("Chain is not contained in any group");
        }

        const auto & worms () const {return _worms;}


        int findWorm(int iChain)
        {
            for ( int i=0;i<_worms.size();i++)
            {
                const auto & wormC = _worms[i];
                if ( (iChain == wormC.iChainHead ) or (iChain == wormC.iChainTail) )
                {
                    return i;
                }
            }

            return -1;
        }

        bool isWorm(int iChain)
        {
            return findWorm(iChain)==-1;
        }

        auto & getMask()  {return _mask;}
        
        protected:

        void updateMask(const worm & wormOld,const worm & wormNew);
        void deleteBeads(   std::array<int,2> timeRange, int iChain ); // deactive the beads in the mask

        void createBeads(   std::array<int,2> timeRange, int iChain ); // activates the bead in the mask

        private:
        int M;
        int N;

        std::vector<particleGroup> particleGroups;
        std::vector<worm> _worms; // order is not preserved during closes
        configurationsStorage_t _data;
        mask _mask;

        int _nParticles;
        particleGroup emptyGroup;
    };

    using configurations_t = pimcConfigurations; 

};

#endif