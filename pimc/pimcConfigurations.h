
#ifndef PIMCCONFIGURATIONS
#define PIMCCONFIGURATIONS

#include <vector>
#include "traits.h"
#include "tools.h"
#include <unsupported/Eigen/CXX11/Tensor>
#include "qmcExceptions.h"
#include "toolsPimc.h"
#include <list>


namespace pimc
{

class maskTensor
{

public:
    maskTensor(){}
    maskTensor( int nChains, int nBeads) : _mask(nChains,nBeads+2) {_mask.setConstant(0);}

    void setConstant(Real value){_mask.setConstant(value);}

    const auto & dataTensor() const {return _mask;}

    const Real * data() const {return _mask.data();}

    Real * data()  {return _mask.data();}
    

    auto & operator()(int i,int j) {return _mask(i,j+1);}
    const auto &  operator()(int i,int j) const {return _mask(i,j+1);}
    void resize(int nChains,int nBeads) {_mask.resize(nChains,nBeads+2);}

private:
    Eigen::Tensor<Real, 2> _mask;
};

enum sector_t{ diagonal = 0 , offDiagonal = 1 , any = 2} ;


struct particleGroup
{
    static const int maxWorms = 1;


    particleGroup( int iStart_,int iEnd_,int iEndExtended_,Real mass_ = 1, sector_t sector_ = sector_t::diagonal) : 
    iStart(iStart_),iEnd(iEnd_),iEndExtended(iEndExtended_),mass(mass_) {
    }
    bool contains(int iParticle) const {return (iParticle>= iStart) and (iParticle<=iEnd);}

    bool contains(const std::array<int,2> & particleRange) const {return contains(particleRange[0]) and contains(particleRange[1]) ;}


    bool isOpen() const {return (heads.size()) > 0 or (tails.size() > 0); }

    
    void pushHead(int iHead) 
    {
        heads.push_back(iHead);
    }
    void pushTail(int iTail) 
    {
        tails.push_back(iTail);
    }
    void popTail(int iTail)
    {
        tails.resize( tails.size() - 1 );
    }
    void popHead(int iHead)
    {
        heads.resize( heads.size() - 1 );
    }

    void removeHead(int iChain)
    {
        auto & _heads = heads;

        auto it = std::find(_heads.begin(),_heads.end(),iChain) ;
        if ( it== _heads.end() )
        {
            throw invalidState("Could not find the head among registered heads.");
        }

    std::swap(*it,*(_heads.end() - 1) );
    _heads.resize(_heads.size() - 1);


    }

    auto range() const {return std::array<int,2> {iStart,iEnd};}

    void removeTail(int iChain)
    {
     
        auto it = std::find(tails.begin(),tails.end(),iChain) ;
        if ( it== tails.end() )
        {
            throw invalidState("Could not find the head among registered heads.");
        }

    std::swap(*it,*(tails.end() - 1) );
    tails.resize(tails.size() - 1);

    }

    bool operator==(particleGroup & groupB)
    {
        return (groupB.iStart == iStart) and (groupB.iEnd == iEnd);
    }

    int iStart; // start of the particle group
    int iEnd; // end of the active group
    int iEndExtended; // extended memory for additional particles
    Real mass; // mass of the particles in the group
    auto size() const {return iEnd - iStart + 1;}
    std::vector<int> heads;
    std::vector<int> tails;
};


class pimcConfigurations;

enum ensamble_t {grandCanonical=2, canonical=1,semiGrandCanonical=3};

struct chain
{

public:
    chain() ;

    bool isOpen() const {return hasHead() or hasTail() ;}
    bool hasHead() const {return next==-1;}
    bool hasTail() const {return prev==-1;}

    bool isEmpty() const {return head -  tail <= 1; }
    void checkChainValid(); 

    int prev; // previous chain , null for a tail
    int next; // next chain, null for a head
    int head; // time Slice of the Head, will be M for closed chains and tails
    int tail; // time slice of the tail , -1 if not a a worm tail
    
};

    class pimcConfigurations
    {
        public:
        using configurationsStorage_t =  Eigen::Tensor<Real, 3> ;
        using tags_t = maskTensor;
        
        void saveHDF5(const std::string & filename);
        static pimcConfigurations loadHDF5(const std::string & filename);

        pimcConfigurations() : pimcConfigurations(0,getDimensions(), {} 
        ) {}

        int nParticles() const ;

        pimcConfigurations(size_t timeSlices, int dimensions, const std::vector<particleGroup> & particleGroups_);
        
        void setChemicalPotential(const std::vector<Real> & mu);
        void setDifferenceChemicalPotential(const std::vector<Real> & deltaMu );


        auto &  dataTensor() {return _data;}
        const auto &  dataTensor() const {return _data;}

        auto data() {return _data.data();}

        const auto data() const {return _data.data();}
        
        auto nChains() const {return N;} 


        auto nBeads() const {return M; } 

        void fillHead(int i);

        void fillHeads();

        const auto & getChain(int i) {return _chains[i];}

        const auto & getChain(int i) const {return _chains[i];}

        void setHead( int iChain, int delta);
        void setTail( int iChain, int delta);
        void setHeadTail( int iChain, int newHead, int newTail);


        void save(const std::string & directoryName,const std::string & format="csv") const;

        void load(const std::string & directoryName);

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
         

          void translateData(const std::array<int,2> & timeRange, const std::array<int,2> & particleRange,const std::array<Real,getDimensions() > & delta);

          void copyData(const std::array<int,2> & timeRange, int iParticleFrom, int timeOffset, int iParticleTo)
         {
             pimcConfigurations::copyData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeOffset,iParticleTo);
         }


        const auto & tails (int iGroup) const {return particleGroups[iGroup].tails;}
        const auto & heads (int iGroup) const {return particleGroups[iGroup].heads;}


        bool isOpen(int iGroup) const { return particleGroups[iGroup].isOpen()  ;}

        bool isOpen() const
        {
            bool open=false;
            for(int i=0;i<particleGroups.size();i++)
            {
                open=open or particleGroups[i].isOpen();
            }
            return open;
        }

        void swapTails(int iChainLeft,int iChainRight);

         void copyDataToBuffer(  Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset=0) const;

         void copyDataFromBuffer(const Eigen::Tensor<Real,2> & buffer, const std::array<int,2> & timeRange, int iParticle ,int timeOffset=0);

        void swapData( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo );

       void swapTags( pimcConfigurations & confFrom, const std::array<int,2> & timeRangeFrom,const std::array<int,2> & particleRangeFrom ,
        pimcConfigurations & confTo,
       int timeOffsetTo, int particleOffsetTo );


       void swapData(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::swapData(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }

         void swapTags(const std::array<int,2> & timeRange, int iParticleFrom, int iParticleTo)
         {
             pimcConfigurations::swapTags(*this,timeRange,{iParticleFrom,iParticleFrom},(*this),timeRange[0],iParticleTo);
         }

        void swapData(int iParticleFrom, int iParticleTo)
         {
             swapData({0,nBeads()},iParticleFrom,iParticleTo);
         }

         void swapTags(int iParticleFrom, int iParticleTo)
         {
             swapTags({0,nBeads()},iParticleFrom,iParticleTo);
         }

        void swap(int iParticleFrom, int iParticleTo);

        void deleteBeads(   std::array<int,2> timeRange, int iChain ); // deactive the beads in the mask

        void createBeads(   std::array<int,2> timeRange, int iChain ); // activates the bead in the mask
        
        const auto & getChemicalPotential() const { return mu;}
        
        const auto & getChemicalPotential(int iSet ) const { return mu[iSet];}


        int nParticles( int i) const;

        const particleGroup & getGroupByChain(int iChain) const {
            for (auto & group : particleGroups )
            {
                if (group.contains(iChain) )
                {
                    return  group;       
                }
            }
             throw invalidInput("Chain " + std::to_string(iChain) + " is not contained in any group");
        }



        void setRandom( const  std::array<Real,DIMENSIONS> & lBox,randomGenerator_t & rand);



        const auto & getTags() const  {return _tags;}

        int pushChain(particleGroup & group);
        int pushChain(int iGroup);


        int next(int iChain);
        int prev(int iChain);
        
        int getHead(int iChain) {return _chains[iChain].head;};
        int getTail(int iChain){return _chains[iChain].tail; };

        void removeChain(int iChain); // orders of chain is not mantained

        void removeChains(int iChain1, int iChain2);

        void join( int iChainLeft, int iChainRight);

        void deleteHeadFromList(int iChain); // deletes the head from the list

        void deleteTailFromList(int iChain); // deletes the tail from the list

        std::list<int> buildPolimerList(int iChain) const; // build a list of all chains in the same permutation cycle as iChain 

        void setEnsamble(ensamble_t ensamble2){ensamble=ensamble2;}
        
        auto getEnsamble() const {return ensamble;}
        
        std::ostream &  operator>>( std::ostream & f) const;

        protected:

         particleGroup & getModifiableGroupByChain(int iChain) {
            for (auto & group : particleGroups )
            {
                if (group.contains(iChain) )
                {
                    return  group;       
                }
            }

            int iEnd=0;

            throw invalidInput("Chain " + std::to_string(iChain) + " is not contained in any group");
        }



        void updateTags(const chain & chainToUpdate);

        
        const auto & getChains() {return _chains; }

        


        int nGroups() const {return particleGroups.size();}

        
        

        private:

        int M;
        int N;


        std::vector<particleGroup> particleGroups;
        configurationsStorage_t _data;
        tags_t _tags;
        std::vector<chain> _chains;

        std::vector<Real> mu;

        ensamble_t ensamble;
        int _nParticles;

    };

    using configurations_t = pimcConfigurations; 


class configurationsSampler
{
    public:
    configurationsSampler() : uniformRealNumber(0,1),normal(0,1),D(0.5) {}
    
    int sampleChain(const configurations_t & confs,randomGenerator_t & randG);

    int sampleChain( const configurations_t & confs,int iGroup,randomGenerator_t & randG);

    int sampleGroup(const configurations_t & confs,randomGenerator_t & randG);




    void sampleFreeParticlePosition(std::array<Real,getDimensions()> & x,const std::array<Real,getDimensions()> & mean,Real tau,randomGenerator_t & randG,Real mass=1);    

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
    
    std::normal_distribution<Real> normal;
    
     Real D ;

};

int nParticlesOnClose(const pimcConfigurations & configurations, int set);


bool checkTimePeriodicBoundaryConditions( const pimcConfigurations & confs, const geometry_t & geo  );

void restrictToBox( pimcConfigurations & confs, const geometry_t & geo);

int nParticlesOnCloseAfterHeadShift(const pimcConfigurations & configurations, int set, int headShift);


struct configurationsRestriction
{
    public:

    virtual bool check(const pimcConfigurations & pimc){return true;};

    private:


};



struct semiCanonicalconfigurationsRestriction : public configurationsRestriction
{
    public:

    semiCanonicalconfigurationsRestriction(std::array<int,2> sets, std::array<int,2> nMin,std::array<int,2> nMax);


    semiCanonicalconfigurationsRestriction( const json_t & j);


    bool check(const pimcConfigurations & confs) override;


    virtual void setHeadShift(int shift, int set);


    std::array<int,2>  nParticlesOnFullClose( const pimcConfigurations & confs);


    private:

    std::array<int , 2> _nMin;
    std::array<int , 2> _nMax;
    std::array<int , 2> _sets;
    std::array<int,2> _shift;

    
};

struct particleRestriction : public configurationsRestriction
{

    particleRestriction(const json_t & j);


    virtual bool check( const pimcConfigurations & confs);

    const auto & nMin() const { return _nMin;}
    const auto & nMax() const { return _nMax;}
    const auto & sets() const { return _sets;}


    protected:
    std::vector<int> _nMin;
    std::vector<int> _nMax;
    std::vector<int> _sets;
    
};


struct wormsOpenRestriction : public particleRestriction
{

    wormsOpenRestriction( const json_t & j ) : particleRestriction::particleRestriction(j)
    {
        auto setToOpen=j["setA"].get<int>();
        for(int i=0;i<_nMin.size();i++)
        {
            if (_sets[i] != setToOpen)
            {
                _nMax[i]-=1;
            }
        }

    }


    virtual bool check( const pimcConfigurations & confs) {return true;}


    private:

};

struct wormsCloseRestriction : public particleRestriction
{

    wormsCloseRestriction( const json_t & j ) : particleRestriction::particleRestriction(j)
    {
        setToClose=j["setA"].get<int>();
        iSetToClose=-1;

        for (int i=0;i<_sets.size();i++)
        {
            if (_sets[i]==setToClose)
            {
                iSetToClose=i;
            }
        }
    }

    virtual bool check( const pimcConfigurations & confs) override;


    private:

    int setToClose;
    int iSetToClose;

};


struct fullCloseRestriction : public particleRestriction
{

    fullCloseRestriction( const json_t & j ) : particleRestriction::particleRestriction(j)
    {
        setToClose=j["setA"].get<int>();
        iSetToClose=-1;

        for (int i=0;i<_sets.size();i++)
        {
            if (_sets[i]==setToClose)
            {
                iSetToClose=i;
            }
        }
    }

    virtual bool check( const pimcConfigurations & confs) override;

    private:

    int setToClose;
    int iSetToClose;

};



struct advanceRestriction : public particleRestriction
{

    advanceRestriction( const json_t & j ) : particleRestriction::particleRestriction(j)
    {
        setA=j["setA"].get<int>();
        iSetA=-1;

        for (int i=0;i<_sets.size();i++)
        {
            if (_sets[i]==setA)
            {
                iSetA=i;
            }
        }
    }

    virtual bool check( const pimcConfigurations & confs) override;

    private:

    int setA;
    int iSetA;

};

struct recedeRestriction : public particleRestriction
{
    recedeRestriction( const json_t & j ) : particleRestriction::particleRestriction(j)
    {
        setA=j["setA"].get<int>();
        iSetA=-1;

        for (int i=0;i<_sets.size();i++)
        {
            if (_sets[i]==setA)
            {
                iSetA=i;
            }
        }


    }


    virtual bool check( const pimcConfigurations & confs) override;

    private:

    int setA;
    int iSetA;

};



};




#endif