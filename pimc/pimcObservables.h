#ifndef PIMC_OBSERVABLES_H
#define PIMC_OBSERVABLES_H

#include "tools.h"
#include "pimcConfigurations.h"
#include "action.h"
#include "accumulators.h"
#include "nConnectedChains.h"
#include<filesystem>

namespace pimc
{

class scalarObservable;
class vectorObservable;

class scalarEstimator
{
    public:
    using observable_t = scalarObservable;

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) = 0;


};

class vectorEstimator
{
    public:
    using observable_t = vectorObservable;
    using accumulator_t = vectorAccumulator<Real>;

    virtual void operator()(configurations_t & configurations, firstOrderAction & S,  accumulator_t & histAcc) = 0;

};

class histogramObservable;

class histogramEstimator
{
    public:
    using accumulator_t = histogramAccumulator<Real>;
    using observable_t = histogramObservable;


    virtual void operator()(configurations_t & configurations, firstOrderAction & S,  accumulator_t & histAcc) = 0;

};

class observable
{
public:

    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)=0;

    virtual void out(size_t iteration)=0;

    virtual void clear()=0;

    virtual bool isValidSector(const configurations_t & confs);      


};

namespace fs = std::filesystem;

class scalarObservable : public observable
{
public:
    using estimator_t = scalarEstimator;

    scalarObservable(std::shared_ptr<scalarEstimator> ob_ , std::string label_,bool append=true) : label(label_),filename(label_ + ".dat"),delim(" "){
        ob=ob_;
        if (append)
        {
            if (fs::exists(fs::path(filename)))
            {
                f.open(filename,std::fstream::app);
            }
            else
            {
                f.open(filename);
                f << "iteration" << delim << label<<std::endl;

            }


        }
        else
        {
            
        }
    }


    scalarObservable(std::shared_ptr<scalarEstimator> ob_ , const json_t & j ) : scalarObservable(ob_,j["label"].get<std::string>() ) {
        
    }

    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)
    {
        acc+=(*ob)(configurations,S);
        
    }

    virtual void out(size_t iteration)
    {
        if ( acc.getWeight() != 0 )
        {
            f << iteration << delim <<  acc.average() << std::endl;
        }

        f << std::flush;
    }



    virtual void clear()
    {
        acc.clear();
    }

    ~scalarObservable()
    {
        f.close();
    }

    Real average()
    {
        return acc.average();
    }

    Real weight()
    {
        return acc.getWeight();
    }



    auto getEstimator() { return ob;}


    private:

    std::shared_ptr<scalarEstimator> ob;
    scalarAccumulator<Real> acc;
    std::ofstream f;
    std::string filename;
    std::string label;
    std::string delim;
};

class histogramObservable : public observable
{
public:

    using estimator_t = histogramEstimator;
    

    histogramObservable(std::shared_ptr<histogramEstimator> ob_ , std::string label_,size_t size,Real min,Real max) : label(label_),filename(label_ + ".dat"),delim(" "){
        ob=ob_;

        if (fs::exists(fs::path(filename)))
        {
            f.open(filename,std::fstream::app);
        }
        else
        {
            f.open(filename);
            f << "iteration" << delim << "x" << delim <<  label<<std::endl;

        }
        
        acc.resize(size,min,max);

    }

     histogramObservable(std::shared_ptr<histogramEstimator> ob_ , const json_t & j) : histogramObservable(ob_,j["label"].get<std::string>() , j["bins"].get<size_t>(),j["minx"].get<Real>(), j["maxx"].get<Real>() ) {}


    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)
    {
        (*ob)(configurations,S,acc);
    }

     auto weight() const
    {
        return acc.weight();
    }

    
    virtual void out(size_t iteration) 
    {
        auto av = acc.average();
        
        if ( weight() != 0 )
        {
            for(int i=0;i<acc.size();i++ )
            {
                f << iteration << delim <<  acc.x(i) << delim << av(i) <<  std::endl;
            }
            
        }
    }


    virtual void clear()
    {
        acc.clear();
    }

    ~histogramObservable()
    {
        f.close();
    }

    auto average() const
    {
        return acc.average();
    }

   
    private:

    histogramAccumulator<Real> acc;
    std::ofstream f;
    std::string filename;
    std::string label;
    std::string delim;
    std::shared_ptr<histogramEstimator> ob;
    
};


class vectorObservable : public observable
{
public:

    using estimator_t = vectorEstimator;

    vectorObservable(std::shared_ptr<vectorEstimator> ob_ , std::string label_ , size_t size ) : label(label_),filename(label_ + ".dat"),delim("\t"){
        ob=ob_;

        if (fs::exists(fs::path(filename)))
            {
                f.open(filename,std::fstream::app);
            }
            else
            {
                f.open(filename);
                f << "iteration" << delim << "x" << delim << label << std::endl;
            }
        
        acc.resize(size,0);
    }

    vectorObservable(std::shared_ptr<vectorEstimator> ob_ , const json_t & j) : vectorObservable( ob_, j["label"].get<std::string>(),j["size"].get<size_t>() ) {}

    virtual void accumulate(configurations_t & configurations, firstOrderAction & S)
    {
        (*ob)(configurations,S,acc);
    }
    
    virtual void out(size_t iteration) 
    {
        auto av = acc.average();

        for(int i=0;i<acc.size();i++ )
        {
                if ( acc.weights()(i) != 0 )
                {
                    f << iteration << delim << i  << delim << av(i) <<  std::endl;
                }
        }

    }


    virtual void clear()
    {
        acc.clear();
    }

    ~vectorObservable()
    {
        f.close();
    }

    auto average() const
    {
        return acc.average();
    }

   
    private:

    vectorAccumulator<Real> acc;
    std::ofstream f;
    std::string filename;
    std::string label;
    std::string delim;
    std::shared_ptr<vectorEstimator> ob;
    
};



class thermodynamicEnergyEstimator : public scalarEstimator 
{

    public:
    
    thermodynamicEnergyEstimator(){}
    
    thermodynamicEnergyEstimator(const json_t & j) : thermodynamicEnergyEstimator() {}


    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);
};

class thermodynamicEnergyEstimatorMagnetization : public vectorEstimator 
{
    public:

    using accumulator_t = vectorAccumulator<Real>;

    thermodynamicEnergyEstimatorMagnetization(){
        setA=0;
        setB=1;
    }

    thermodynamicEnergyEstimatorMagnetization(const json_t & j) : energyEst(j)   {setA=0;setB=1;}

    virtual void operator()(configurations_t & configurations, firstOrderAction & S, accumulator_t & acc);

    private:
    int setA;
    int setB;
    thermodynamicEnergyEstimator energyEst;
};


class superfluidFractionEstimator : public vectorEstimator
{
    public:

    using accumulator_t = vectorAccumulator<Real>;

    superfluidFractionEstimator(){};

    superfluidFractionEstimator(const json_t & j) {};

    virtual void operator()(configurations_t & configurations, firstOrderAction & S, accumulator_t & acc);

    private:

};

class particleNumberEstimator : public scalarEstimator 
{
    
    public:
    
    particleNumberEstimator( const std::vector<int> & sets_) : _sets(sets_) {}

    particleNumberEstimator(int setA) : _sets({setA}) {}


    
    particleNumberEstimator(const json_t & j) : particleNumberEstimator( j["sets"].get<std::vector<int> >() ) {}
    

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) {
        int N=0;
        for( auto set : _sets)
        {
            N+=configurations.nParticles(set);
        }
    return N;
    }

    private: 
    std::vector<int> _sets;
};


class particleNumberSquaredEstimator : public scalarEstimator 
{
    
    public:
    
    particleNumberSquaredEstimator( const std::vector<int> & sets_) : _sets(sets_) {}

    particleNumberSquaredEstimator(int setA) : _sets({setA}) {}

    
    particleNumberSquaredEstimator(const json_t & j) : particleNumberSquaredEstimator( j["sets"].get<std::vector<int> >() ) {}
    
        

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) {
        int N=0;
        for( auto set : _sets)
        {
            N+=configurations.nParticles(set);
        }
    return N*N;
    }

    private: 
    std::vector<int> _sets;
};




class magnetizationSquaredEstimator : public scalarEstimator 
{
    
    public:
    
    magnetizationSquaredEstimator(int setA,int setB) : _setA(setA), _setB(setB) {}
    
    magnetizationSquaredEstimator(const json_t & j) :  magnetizationSquaredEstimator( j["groupA"].get<int>() , j["groupB"].get<int>()  ) {}
    
    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) {

        int NA = configurations.nParticles(_setA) ;
        int NB = configurations.nParticles(_setB) ;
        
        return (NA - NB)*(NA-NB);
        }
        
    private: 
    
    int _setA;
    int _setB;

};

class magnetizationEstimator : public scalarEstimator 
{
    
    public:
    
    magnetizationEstimator(int setA,int setB) : _setA(setA), _setB(setB) {}
    
    magnetizationEstimator(const json_t & j) :  magnetizationEstimator( j["groupA"].get<int>() , j["groupB"].get<int>()  ) {}
    

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) {

        int NA = configurations.nParticles(_setA) ;
        int NB = configurations.nParticles(_setB) ;

        return std::abs(NA - NB);
        }
    
    
    private: 

    int _setA;
    int _setB;

};




class lengthEstimator : public scalarEstimator 
{
    
    public:

    lengthEstimator( );

    lengthEstimator( int l);
    

    lengthEstimator(const json_t & j) :  lengthEstimator( j["chain"].get<int>()   ) { }
    
    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);

    void setStartFromHeadOrTail(int iGroup);

    void setStartingBead(int iChain,int t0) {_iChain=iChain;_startFromHeadOrTail=false;_t0=t0;}


    private:

    int  initialChain( const configurations_t & configurations) const ;

    int initialTime( const configurations_t & configurations) const ;

    private: 

    int _iChain;
    int _t0;
    int _l;
    bool _startFromHeadOrTail;
    int _iGroup;
    bool _setMaxLength;

};


class closedLengthEstimator : public scalarEstimator 
{
    
    public:

    closedLengthEstimator( int set ) : _set(set) {}

    closedLengthEstimator(const json_t & j) :  closedLengthEstimator( j["set"].get<int>()   ) { }

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);

    private: 

    int _set;

};

class nBeadsInWormEstimator : public scalarEstimator 
{
    public:
    
    nBeadsInWormEstimator( int iGroup ) : _iGroup(iGroup) {};
    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);


    private:
    int _iGroup;

};

class nConnectedChainsEstimator : public scalarEstimator
{
    public:

    nConnectedChainsEstimator( int set_) : _set(set_){};
    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);


    private:

    nConnectedChains chainCounter;
    int _set;
};

class virialEnergyEstimator : public scalarEstimator
{
    public:
  
    virialEnergyEstimator(int nMax, int MMax) : buffer(nMax,getDimensions(),MMax+1),rC(nMax,getDimensions(),MMax) {}


    virialEnergyEstimator(const json_t & j) : virialEnergyEstimator(j["nChains"].get<int>(),j["nBeads"].get<int>() ) {}

        
    
    Real operator()(configurations_t & configurations, firstOrderAction & S);
    
    private:
    Eigen::Tensor<Real,3> buffer;
    Eigen::Tensor<Real,3> rC;
};

class pairCorrelation : public histogramEstimator
{
    public:

    pairCorrelation(int setA,int setB);

    pairCorrelation(const json_t & j) : pairCorrelation(j["setA"
    ].get<int>() , j["setB"].get<int>()  ) {}




    void operator()(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc);

    private:

    Real getNormalizationFactor(const configurations_t & configurations, const firstOrderAction & S , const accumulator_t & acc) const ;

    void accumulateDistinguishable(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc);

    void accumulateUnDistinguishable(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc);


    int setA;
    int setB;
    std::vector<double> buffer;

};





class oneBodyEstimator : public histogramEstimator
{
    public:

    oneBodyEstimator(int set) : _set(set) {};

    oneBodyEstimator(const json_t & j) : oneBodyEstimator( j["set"].get<int>()  ) {}

    void operator()(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc);

    private:

    int _set;
    
};

class oneBodyObservable : public histogramObservable
{
public:

    using estimator_t = histogramEstimator;

    virtual bool isValidSector(const configurations_t & confs) override;


    oneBodyObservable(std::shared_ptr<oneBodyEstimator> ob_ , std::string label_,size_t size,Real min,Real max,int set) : histogramObservable(ob_,label_,size,min,max) { _set=set;};
    
    oneBodyObservable( const json_t & j) :
    histogramObservable( std::make_shared<oneBodyEstimator>(j["set"].get<int>()) , j)
     {
       _set=j["set"].get<int>(); 
    };




    
   
    private:

    int _set;


};



class angleEstimator : public histogramEstimator
{
    public:

    angleEstimator(int setA,int setB);

    angleEstimator(const json_t & j) : angleEstimator(j["setA"
    ].get<int>() , j["setB"].get<int>()  ) {}




    void operator()(configurations_t & configurations, firstOrderAction & S,accumulator_t & acc);

    private:
    
    int setA;
    int setB;
    std::vector<double> buffer;

};



class magnetizationDistribution : public observable
{
    public:

    magnetizationDistribution( const json_t & j);
    
    virtual void accumulate(configurations_t & configurations, firstOrderAction & S) override;

    virtual void out(size_t t) override;

    void clear() override;

    ~magnetizationDistribution();

    private:
    
    std::vector<Real> _Ms;
    size_t n;
    int _mMin;
    int _mMax;
    int setA;
    int setB;
    std::ofstream f;
    std::string _label;
    bool recordAbsoluteValue;  
    
};

class openRatio 
{
    public:
    openRatio( int nSets );

    void accumulate(configurations_t & configurations);

    auto average( int i) { return nOpen[i]/n  ;  }

    void clear() { n=0; std::fill(nOpen.begin(), nOpen.end(), 0); }

    void out( int i );


    ~openRatio() ;

    private:

    std::vector< Real > nOpen;
    Real n;

    std::string filename;
    std::ofstream f;

};

class virialEnergyEstimatorMagnetization : public vectorEstimator 
{
    public:

    using accumulator_t = vectorAccumulator<Real>;


    virialEnergyEstimatorMagnetization(const json_t & j) : energyEst(j)   {setA=0;setB=1;}

    virtual void operator()(configurations_t & configurations, firstOrderAction & S, accumulator_t & acc);

    private:
    int setA;
    int setB;
    virialEnergyEstimator energyEst;
};




};


#endif