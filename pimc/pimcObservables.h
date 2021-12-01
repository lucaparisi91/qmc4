#ifndef PIMC_OBSERVABLES_H
#define PIMC_OBSERVABLES_H


#include "tools.h"
#include "pimcConfigurations.h"
#include "action.h"
#include "accumulators.h"


namespace pimc
{

class scalarObservable;

class scalarEstimator
{
    public:
    using observable_t = scalarObservable;

    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) = 0;

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

};


class scalarObservable : public observable
{
public:
    using estimator_t = scalarEstimator;
    
    scalarObservable(std::shared_ptr<scalarEstimator> ob_ , std::string label_) : label(label_),filename(label_ + ".dat"),delim(" "){
        ob=ob_;
        f.open(filename,std::fstream::app);
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
        f.open(filename,std::fstream::app);
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
                f << iteration << delim <<  acc.x(i) << delim << av(i)/( acc.x(i)*acc.x(i) ) <<  std::endl;
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



class thermodynamicEnergyEstimator : public scalarEstimator 
{
    
    public:
    
    thermodynamicEnergyEstimator(){}
    
    thermodynamicEnergyEstimator(const json_t & j) : thermodynamicEnergyEstimator() {}


    virtual Real operator()(configurations_t & configurations, firstOrderAction & S);
};


class particleNumberEstimator : public scalarEstimator 
{
    
    public:
    
    particleNumberEstimator(int setA) : _setA(setA) {}
    
    particleNumberEstimator(const json_t & j) : particleNumberEstimator( j["set"].get<int>() ) {}


    virtual Real operator()(configurations_t & configurations, firstOrderAction & S) {return configurations.nParticles(_setA);}

    private: 
    int _setA;
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
        
        return std::abs(NA - NB)*1./(NA+NB);
        }

    private: 

    int _setA;
    int _setB;

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


}


#endif