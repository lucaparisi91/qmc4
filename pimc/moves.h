#ifndef MOVESPIMC_H
#define MOVESPIMC_H

#include "traits.h"
#include "pimcConfigurations.h"
#include "tools.h"
#include "metropolis.h"
#include "geometryPMC.h"
#include "towerSampler.h"
#include "toolsPimc.h"
#include <cassert>
#include "levyReconstructor.h"


namespace pimc
{

    class firstOrderAction;
    class action;

class timeSliceGenerator
{
    public:
    timeSliceGenerator(){}

    std::array<int, 2> operator()(randomGenerator_t & randG, int nBeads, int maxBeadLength);

    private:
    std::uniform_real_distribution<float> uniformRealNumber;
}; 

class move 
{
    public:

    move(const std::vector<int> & sets ) : _sets(sets){}

    virtual bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)=0;

    const auto & getSets() const {return _sets; }

    private:
    std::vector<int> _sets;
};

class singleSetMove : public move
{
    public:
    singleSetMove(int set) : move({set}){}
    singleSetMove(const json_t & j) : singleSetMove::singleSetMove(j["group"].get<int>() ){}

    virtual bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)=0;

    int getSet(){return getSets()[0];}

    auto & getChainSampler() {return _chainSampler;}

    private:
    configurationsSampler _chainSampler;
};

class twoSetMove : public move
{
    public:
    twoSetMove(int setA, int setB) : move({setA,setB}){}
    twoSetMove(const json_t & j) : twoSetMove::twoSetMove(j["setA"].get<int>() ,j["setB"].get<int>() ){}
    

    virtual bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG)=0;


    int getSetA(){return getSets()[0];}

    int getSetB(){return getSets()[1];}

    auto & getChainSampler() {return _chainSampler;}

    private:
    configurationsSampler _chainSampler;
};





class sectorTableMoves
{
    public:
    
    sectorTableMoves() : totalWeight(0){}
    
    void push_back( move * move_,Real weight,const std::string & name="Unkown move") ; 

    bool attemptMove(configurations_t & confs,firstOrderAction & S,randomGenerator_t & randG);

    Real acceptanceRatio(int i) {return _nSuccess[i]*1./_nTrials[i];}

    std::ostream & operator>> (std::ostream & os);

    void resetCounters();

    int size() const {return _moves.size();}

    const auto &  nTrials() const  {return _nTrials;}

    const auto &  nSuccess() const  {return _nSuccess;}
    
    const auto &  names() const  {return _names;}

    const auto & getTotalWeight() const {return totalWeight;}

    const auto & getNames() const {return _names;}

    const auto & getWeights() const {return _weights;}


    private:



    int sample(randomGenerator_t & random); 
    std::vector<move*> _moves;

    std::vector<Real> _nTrials;
    std::vector<Real> _nSuccess;
    std::vector<std::string> _names;

    std::vector<Real>  accumulatedWeights;
    std::vector<Real>  _weights;
    
    Real totalWeight;
    towerSampler sampler;
};

std::string reportMoves(const sectorTableMoves & tab);


class tableMoves
{
    public:

    tableMoves(){}

    void push_back( move * move_,Real weight,sector_t sector,const std::string & name="Unkown move") ; 
    

    bool attemptMove(configurations_t & confs, firstOrderAction & S,randomGenerator_t & randG);


    std::ostream & operator>> (std::ostream & os);

    void resetCounters();
    
    auto & getOpenSectorTable(int i) {return openTabs[i];}
    const auto & getClosedSectorTable(int i) const {return closedTabs[i];}


    private:

    void resize(int nSets);

    int nRegisteredSets() const {return openTabs.size();}


    std::vector<sectorTableMoves> openTabs;
    std::vector<sectorTableMoves> closedTabs;

    std::vector<Real> nOpenSectorMoves;
    std::vector<Real> nClosedSectorMoves;
    configurationsSampler _chainSampler;
};


class levyMove : public singleSetMove
{
    public:

    levyMove( int maxBeadLength, int set);
    levyMove(const json_t & j) : levyMove::levyMove(j["reconstructionMaxLength"].get<int>() , j["set"].get<int>()  )
    {}



    bool attemptMove(configurations_t & confs , firstOrderAction & S, randomGenerator_t & randG);

    private:

    bool isValidSlice( configurations_t & confs , const std::array<int,2> timeRange,int iChain) const;

    int maxBeadLength;

    levyReconstructor _levy;
    std::uniform_real_distribution<float> uniformRealNumber;
    
    metropolis sampler;
    configurationsSampler confsSampler;
    Eigen::Tensor<Real,2> buffer;
    timeSliceGenerator tGen;

};


class openMove : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    openMove(Real C_ , int set,int maxReconstructedLength_=1) ;

    openMove(const json_t & j) : openMove(j["C"].get<Real>() ,j["set"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) {}


    void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0);};
    void setLengthCut( int l1){setLengthCutRandom=false;lengthCut=l1;}

    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    
    bool attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    
    bool attemptSemiGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    

    bool attemptCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    private:

    Real openCloseRatioCoefficient(int N,int M);


    Real C;
    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;


    int startingBead;
    bool setStartingBeadRandom;
    bool setLengthCutRandom;
    bool setStartingChainRandom;
    int lengthCut;
    int startingChain;
};


class semiOpenMove : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    semiOpenMove(Real C , int setA ,int maxLength_=1 ) ;

    semiOpenMove(const json_t & j) : semiOpenMove(j["C"].get<Real>() ,j["group"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) {}

    void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0);};
    void setLength( int l1){setLengthCutRandom=false;lengthCut=l1;}

    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    
    private:

    Real openCloseRatioCoefficient(int N,int M);

    Real C;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int length;

    Eigen::Tensor<Real,2> buffer;


    int startingBead;
    int startingChain;
    bool setStartingBeadRandom;
    bool setLengthCutRandom;
    bool setStartingChainRandom;
    int lengthCut;
    int _maxLength;
};

class fullSemiCanonicalOpenMove : public twoSetMove
{
    public:

    fullSemiCanonicalOpenMove(Real CA , int setA , int setB, int maxLength_=1 );

    fullSemiCanonicalOpenMove( const json_t & j) : fullSemiCanonicalOpenMove(j["C"].get<Real>() ,j["setA"].get<int>() ,j["setB"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) {}


    void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0);};

    void setLength( int l1){setLengthCutRandom=false;lengthCut=l1;}

    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    
    private:

    Real openCloseRatioCoefficient(int N,int M);

    Real C;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int length;

    Eigen::Tensor<Real,2> buffer;


    int startingBead;
    int startingChain;
    bool setStartingBeadRandom;
    bool setLengthCutRandom;
    bool setStartingChainRandom;
    int lengthCut;
    int _maxLength;
};

class fullSemiCanonicalCloseMove : public twoSetMove
{
    public:
    
    fullSemiCanonicalCloseMove(Real CA , int setA , int setB, int maxLength_=1 );

    fullSemiCanonicalCloseMove( const json_t & j) : fullSemiCanonicalCloseMove(j["C"].get<Real>() ,j["setA"].get<int>() ,j["setB"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) {}


    void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0);};

    void setLength( int l1){setLengthCutRandom=false;lengthCut=l1;}

    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    
    private:

    Real openCloseRatioCoefficient(int N,int M);

    Real C;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int length;

    Eigen::Tensor<Real,2> buffer;

    int startingBead;
    int startingChain;
    bool setStartingBeadRandom;
    bool setLengthCutRandom;
    bool setStartingChainRandom;
    int lengthCut;
    int _maxLength;
};



class semiCloseMove : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    semiCloseMove(Real C , int setA ,int maxLength_=1 ) ;

    semiCloseMove(const json_t & j) : semiCloseMove(j["C"].get<Real>() ,j["group"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) {}

    void setLength( int l1){setLengthRandom=false;length=l1;}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

     void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0);};


    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};
    
    
    private:

    Real openCloseRatioCoefficient(int N,int M);

    Real C;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int length;

    Eigen::Tensor<Real,2> buffer;


    int startingBead;
    int startingChain;
    bool setStartingBeadRandom;
    bool setLengthRandom;
    bool setStartingChainRandom;
    int _maxLength;

};


class openMoveTest : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    openMoveTest(Real C_ , int set,int maxLength=1, int startingBead=0) ;

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setOrder(int order_){order=order_;}

    private:


    Real C;
    int _maxLength;
    int order;


    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int startingBead;

};

class closeMoveTest : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    closeMoveTest(Real C_ , int set,int maxLength=1, int startingBead=0) ;

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
     void setOrder(int order_){order=order_;}
    private:


    Real C;
    int _maxLength;
    int order;


    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    int startingBead;

};


class createWorm : public singleSetMove
{
    public:
    // splits a chain in two morms with one overlapping bead
    createWorm(Real C_ , int set,int maxReconstructedLength_=1,Real sigma_=1);


    createWorm(const json_t & j) : createWorm(j["C"].get<Real>() ,j["set"].get<int>() ,j["reconstructionMaxLength"].get<int>() ) 
    {
        if ( j.find("alpha") != j.end() )
        {
            Real sigma = 1/sqrt(2*j["alpha"].get<Real>() );
            setSigma(sigma);
        }

        if ( j.find("firstParticleDistribution") != j.end() )
        {
            std::string dis = j["firstParticleDistribution"].get<std::string>();

            if (dis == "uniform")
            {
                uniformDistributionFirstParticle=true;
            }
            else if (dis == "gaussian")
            {
                uniformDistributionFirstParticle=false;
            }
            else
            {
                throw std::runtime_error("Unkown first particle distribution: " + dis);

            }
            


        }


    }


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    Real probabilityInitialPosition(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
    {
        if (not uniformDistributionFirstParticle)
        {
            return probabilityInitialPositionGaussian(geo,x);
        }
        else
        {
            return probabilityInitialPositionUniform(geo);
        }
    }


    void sampleFirstParticle(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG)
    {
        if (uniformDistributionFirstParticle)
        {
            sampleFirstParticleUniform(geo,x1,randG);
        }
        else
        {
            sampleFirstParticleGaussian(geo,x1,randG);
        }

    }


    void sampleFirstParticleGaussian(const geometry_t & geo , std::array<Real,getDimensions()> & x1 , randomGenerator_t & randG)
    {
        for(int d=0;d<getDimensions();d++)
        {
            x1[d]=gauss(randG)*std::sqrt(sigma2);
        }

    }


    void sampleFirstParticleUniform(const geometry_t & geo , std::array<Real,getDimensions()> & x1 ,randomGenerator_t & randG)
    {
        for(int d=0;d<getDimensions();d++)
        {
            x1[d]=(uniformRealNumber(randG)-0.5)*geo.getLBox(d);
        }

    }

    Real probabilityInitialPositionGaussian(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
    {
        Real sum=0;

        for(int d=0;d<getDimensions();d++)
        {
            sum+= -0.5*log(2*M_PI*sigma2) -0.5 * x[d]*x[d] / sigma2;    
        }
        return sum;
    }

    Real probabilityInitialPositionUniform(const geometry_t & geo) const
    {
        Real V=1;
        for(int d=0;d<getDimensions();d++)
        {
            V*=geo.getLBox(d)  ;  
        }
        return -log(V);
    }


    int sampleInitialBead(randomGenerator_t & randG, int M)  ;
    int sampleLength(randomGenerator_t & randG)  ;

    Real getCoefficientAccept(Real C_,int M_,int l_) const
    {
        Real coeff = randomInitialBead ? C_*M_ : C_;
        coeff *= randomLength ? l_  : 1;
        return coeff;
    }


    void setStartingBead(int t0)
    {
        randomInitialBead=false;
        initialBead=t0;
    }


    void setFixedLength()
    {
        randomLength=false;
    }

    void setSigma(Real sigma)
    {
        sigma2=sigma*sigma;
    }

    void setUniformSampling() {uniformDistributionFirstParticle=true;}

    void setGaussianSampling() {uniformDistributionFirstParticle=false;}



    private:
    Real C;
    int _maxReconstructedLength;
    std::array<Real, 3> tmp;
    Real sigma2;

    bool randomLength;
    bool randomInitialBead;
    int initialBead;


    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool uniformDistributionFirstParticle;


};

class deleteWorm : public singleSetMove
{
    public:
    // splits a chain in two worms with one overlapping bead
    deleteWorm(Real C_ , int set,int maxReconstructedLength_=1,Real sigma_=1);


    deleteWorm(const json_t & j) : deleteWorm(j["C"].get<Real>() ,j["set"].get<int>() ,j["reconstructionMaxLength"].get<int>()  )
    {
        
        if ( j.find("alpha") != j.end() )
        {
            Real sigma = 1/sqrt(2*j["alpha"].get<Real>() );
            setSigma(sigma);
        }

        if ( j.find("firstParticleDistribution") != j.end() )
        {
            std::string dis = j["firstParticleDistribution"].get<std::string>();

            if (dis == "uniform")
            {
                uniformDistributionFirstParticle=true;
            }
            else if (dis == "gaussian")
            {
                uniformDistributionFirstParticle=false;
            }
            else
            {
                throw std::runtime_error("Unkown first particle distribution: " + dis);

            }
            


        }


         

    }


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    Real probabilityInitialPosition(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
    {
        if (not uniformDistributionFirstParticle)
        {
            return probabilityInitialPositionGaussian(geo,x);
        }
        else
        {
            return probabilityInitialPositionUniform(geo);
        }
    }


    Real probabilityInitialPositionGaussian(const geometry_t & geo,const std::array<Real,getDimensions()> & x) const
    {
        Real sum=0;

        for(int d=0;d<getDimensions();d++)
        {
            sum+= -0.5*log(2*M_PI*sigma2) -0.5 * x[d]*x[d] / sigma2;    
        }
        return sum;
    }

    Real probabilityInitialPositionUniform(const geometry_t & geo) const
    {
        Real V=1;
        for(int d=0;d<getDimensions();d++)
        {
            V*=geo.getLBox(d)  ;  
        }
        return -log(V);
    }



    Real getCoefficientAccept(Real C_,int M_,int l_) const
    {
        Real coeff = randomInitialBead ? C_*M_ : C_;
        coeff *= randomLength ? l_  : 1;
        return coeff;
    }
    
    void setStartingBead(int t0)
    {
        randomInitialBead=false;
    }

    void setFixedLength()
    {
        randomLength=false;
    }

    void setSigma(Real sigma)
    {
        sigma2=sigma*sigma;
    }

      void setUniformSampling() {uniformDistributionFirstParticle=true;}

    void setGaussianSampling() {uniformDistributionFirstParticle=false;}


    private:
    Real C;
    int _maxReconstructedLength;
    std::array<Real, 3> tmp;
    Real sigma2;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    metropolis sampler;
    bool randomLength;
    bool randomInitialBead;
    bool uniformDistributionFirstParticle;

};


class closeMove : public singleSetMove
{
    public:
    
    // splits a chain in two morms with one overlapping bead
    closeMove(Real C_ , int set,int maxReconstructionLength=1) ;
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    bool attemptCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    bool attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    bool attemptSemiGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    



    closeMove(const json_t & j) : closeMove(j["C"].get<Real>() , j["set"].get<int>(),j["reconstructionMaxLength"].get<int>()   ) {}

    void setStartingBead(int m){setStartingBeadRandom=false; startingBead=m;assert(m>=0 );};
    void setLengthCut( int l1){setLengthCutRandom=false;lengthCut=l1;}

    
    void setStartingChain(int m){setStartingChainRandom=false; startingChain=m;assert(m>=0);};

    void setAtomNumberOnClose(size_t N,const std::vector<int> & sets);




    private:

    bool checkConstraintsOnClose(const configurations_t & confs) ;

    Real openCloseRatioCoefficient(int N,int M);

    Real C;
    std::array<Real, 3> tmp;
    int _maxLength;

    levyReconstructor _levy;
    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;

    int startingBead;
    bool setStartingBeadRandom;
    bool setLengthCutRandom;
    int lengthCut;
    bool setStartingChainRandom;
    int startingChain;

    bool enforceAtomeNumberOnClose;
    size_t atomNumberOnClose;
    std::vector<int> setsMaxParticleNumber;



};

// advance and recede in opposite directions
class moveHead : public singleSetMove
{
    public:
    moveHead(int maxAdvanceLength_,int set);

    moveHead(const json_t & j) : moveHead(j["reconstructionMaxLength"].get<int>() ,j["set"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    private:

    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};

class advanceHead : public singleSetMove
{
    public:
    advanceHead(int maxAdvanceLength_,int set);


    advanceHead(const json_t & j) : advanceHead(j["reconstructionMaxLength"].get<int>() ,j["set"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

     bool attemptMoveSemiCanonical(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    void setFixedLength() {_setRandomLength=false;}

    void setRandomLength() {_setRandomLength=true;}

    void setMaximumParticleNumber(int nMax) {_nMax=nMax;enforceMaxParticleNumber=true;}


    private:

    int sampleLength(randomGenerator_t & randG);


    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool _setRandomLength;
    int _nMax;
    bool enforceMaxParticleNumber;
    
};

class advanceHeadTail : public twoSetMove
{
    public:

    advanceHeadTail(int maxAdvanceLength_,int setA, int setB );

    advanceHeadTail(const json_t & j) : advanceHeadTail::advanceHeadTail(j["reconstructionMaxLength"].get<int>() ,j["setA"].get<int>(), j["setB"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setRandomLength() {_setRandomLength=true;}
    void setFixedLength() {_setRandomLength=false;}

    void setMaximumParticleNumber(int nMax) {_nMax=nMax;enforceMaxParticleNumber=true;}

    private:

    int sampleLength(randomGenerator_t & randG);


    int _maxAdvanceLength;
    std::array<Real, 3> tmp;


    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool _setRandomLength;

    int _nMax;
    bool enforceMaxParticleNumber;

};



class recedeHeadTail : public twoSetMove
{
    public:

    recedeHeadTail(int maxAdvanceLength_,int setA, int setB );

    recedeHeadTail(const json_t & j) : recedeHeadTail::recedeHeadTail(j["reconstructionMaxLength"].get<int>() ,j["setA"].get<int>(), j["setB"].get<int>() ) {}


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    void setRandomLength() {_setRandomLength=true;}
    void setFixedLength() {_setRandomLength=false;}
    
    void setMaximumParticleNumber(int nMax) {_nMax=nMax;enforceMaxParticleNumber=true;}


    private:

    int sampleLength(randomGenerator_t & randG);


    int _maxAdvanceLength;
    std::array<Real, 3> tmp;


    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool _setRandomLength;

    int _nMax;
    bool enforceMaxParticleNumber;

};



class recedeTail : public singleSetMove
{
    public:
    recedeTail(int maxAdvanceLength_,int set);

    recedeTail(const json_t & j) : recedeTail(j["reconstructionMaxLength"].get<int>() ,j["set"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setFixedLength() {_setRandomLength=false;}

    void setRandomLength() {_setRandomLength=true;}


    private:


    int sampleLength(randomGenerator_t & randG);


    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool _setRandomLength;
};

class advanceTail : public singleSetMove
{
    public:
    advanceTail(int maxAdvanceLength_,int set);

    advanceTail(const json_t & j) : advanceTail(j["reconstructionMaxLength"].get<int>() ,j["set"].get<int>() ) {}
    
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setFixedLength() {_setRandomLength=false;}

    void setRandomLength() {_setRandomLength=true;}


    private:


    int sampleLength(randomGenerator_t & randG);


    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    bool _setRandomLength;
};



class advanceHeadTest : public singleSetMove
{
    public:
    advanceHeadTest( int l);


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

   

    private:

    int l;
    const Real D = 0.5;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    configurationsSampler confsSampler;

    metropolis sampler;
};


class recedeHeadTest : public singleSetMove
{
    public:
    recedeHeadTest( int l);

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    private:
    
    int l;
    const Real D = 0.5;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
};

class recedeHead : public singleSetMove
{
    public:
    recedeHead(int maxAdvanceLength_,int set);

    void setFixedLength() {_setRandomLength=false;}

    void setRandomLength() {_setRandomLength=true;}
    

    recedeHead(const json_t & j) : recedeHead(j["reconstructionMaxLength"].get<int>() ,j["set"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setMinParticleNumber( const int & new_nMin_) { _nMin=new_nMin_ ;enforceMinParticleNumber=true;};  
    

    private:
    int sampleLength(randomGenerator_t & randG);

    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;

    bool _setRandomLength;
    int _nMin;
    bool enforceMinParticleNumber;

};

class moveTail : public singleSetMove
{
    public:
    moveTail(int maxAdvanceLength_, int set);

    moveTail(const json_t & j) : moveTail(j["reconstructionMaxLength"].get<int>() , j["set"].get<int>() ) {}

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);
    
    private:


    int _maxReconstructedLength;
    std::array<Real, 3> tmp;

    const Real D = 0.5;
    configurationsSampler confsSampler;
    std::normal_distribution<Real> gauss;
    std::uniform_real_distribution<float> uniformRealNumber;
    levyReconstructor _levy;
    metropolis sampler;
    Eigen::Tensor<Real,2> buffer;
};

class sampleSingleClosedHarmonicChain : public singleSetMove
{
    public:
    sampleSingleClosedHarmonicChain(Real var,int iChain,int M);

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);


    private:
    Real var;
    levyReconstructor _levy;
    metropolis sampler;
    std::normal_distribution<Real> gauss;
    int iChain;
    Eigen::Tensor<Real, 3 > buffer;

};


class sampleSingleOpenHarmonicChain : public singleSetMove
{
    public:
    sampleSingleOpenHarmonicChain(Real var,int iChain,int M);

    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    private:
    Real var;
    levyReconstructor _levy;
    metropolis sampler;
    std::normal_distribution<Real> gauss;
    int iChain;
    Eigen::Tensor<Real, 3 > buffer;

};





// advance and recede in opposite directions
class translateMove : public singleSetMove
{
    public:
    translateMove(Real max_delta, int maxBeads, int set);


    translateMove(const json_t & j) : translateMove(j["delta"].get<Real>() ,(j["nBeads"].get<int>() + 1 ) * getTotParticlesFromJson(j)  , j["set"].get<int>()  ) {}
    
    
    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    private:

    Real _max_delta;
    std::uniform_real_distribution<Real> distr;
    Eigen::Tensor<Real,2> buffer;
    std::array<Real,getDimensions()> delta;
    configurationsSampler confSampler;
    metropolis sampler;
    std::vector<int> currentPolimerList;
};

class swapMove : public singleSetMove
{
    public:
    using geometry_t = geometryPBC_PIMC;
    
    swapMove( int maxStepLength_, int maxN, int set);

    swapMove( const json_t & j ) : swapMove(
        j["reconstructionMaxLength"].get<int>() ,
          j["particles"][ j["set"].get<int>()  ].get<int>() ,  j["set"].get<int>()     ) {}


    bool attemptMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG); // will attempt to perform a new move

    bool attemptGrandCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG); // will attempt to perform a new move

    bool attemptTestMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG); // will attempt to perform a new move
    

    bool attemptCanonicalMove(configurations_t & confs , firstOrderAction & S,randomGenerator_t & randG);

    void setRandomLength()
    {
        _setRandomLength=true;
    }

    void setFixedLength()
    {
        _setRandomLength=false;
    }

    int sampleLength()
    {
        if (_setRandomLength)
        {
            return std::floor(uniformRealNumber(randG) * maxStepLength) + 1;
        }
        else
        {
            return maxStepLength;
        }
       
    }


    private:

    int maxStepLength;
    levyReconstructor _levy;
    randomGenerator_t randG;
    std::uniform_real_distribution<float> uniformRealNumber;
    Eigen::Tensor<Real,2> buffer;
    std::vector<Real> particleSelectionAccWeights;
    Real particleSelectionWeight;
    const Real D = 0.5;
    metropolis metropolisSampler;  
    timeSliceGenerator tGen;
    configurationsSampler confSampler;
    towerSampler particleSampler;
    bool _setRandomLength;
};

template<class T>
move* __createMove(const json_t & j) {return new T(j);}

class moveConstructor
{
    public:

    moveConstructor(std::vector<int> nMaxParticles, int nBeadsMax) : _nMaxParticles(nMaxParticles),_nBeadsMax(nBeadsMax) {}
    
    move* createMove(const json_t & jOuter)
    {
        // creates a copy of the json input and supplement with additional essential information
        json_t j(jOuter);
        j["particles"]=_nMaxParticles;
        j["nBeads"]=_nBeadsMax;


        std::string key=j["kind"].get<std::string>(); 
        creatorMap_t::const_iterator i;
        i=creatorMap.find(key);
        if (i!= creatorMap.end())
        {
	    return (i->second)( j  );
        }
         else
        {
	    throw factoryIdNotRecorded(key);
        }
    }

    template<class T>
    void registerMove(const std::string & key)
    {
        creatorMap[key]= &__createMove<T> ; 
    }

    auto createTable(const json_t & jTable)
    {
        tableMoves tab;

        for (auto & jMove : jTable )
        {
            std::vector<std::string > sectors = jMove["sectors"];
            Real weight = jMove["weight"].get<Real>();

            const auto sets = jMove["sets"].get<std::vector<int> >();

            for ( auto set : sets)
            {

                json_t moveJInput (jMove["move"]);
                moveJInput["set"]=set;


                auto move = createMove(moveJInput);

                std::string kind = jMove["move"]["kind"].get<std::string>();
                
                for (auto sector : sectors)
                {
                    sector_t currentSector;

                    if (sector == "open" )
                    {
                        currentSector=sector_t::offDiagonal;
                    }
                    else if ( sector == "closed")
                    {
                        currentSector = sector_t::diagonal;
                    }
                    else
                    {
                        throw invalidInput("Unkown sector type");
                    } 

                    tab.push_back(move,weight,currentSector,kind);
                }

            }
                 

                
            }
            

        

        return tab;

    }


    private:
    std::vector<int> _nMaxParticles;
    int _nBeadsMax;


    typedef move* (*moveCreatorFunc) ( const json_t & j);
using creatorMap_t = std::map<std::string,moveCreatorFunc>;
    creatorMap_t creatorMap;
};

//moveConstructor::creatorMap =  moveConstructor::creatorMap_t{} ;



}





#endif