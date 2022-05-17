#ifndef ACTION_TWO_BODY_H
#define ACTION_TWO_BODY_H

#include "memory"
#include "kernels.h"
#include "action.h"

namespace pimc
{
//using configurations_t = pimcConfigurations;


struct twoBodyEvaluationPolicy
{

     virtual Real evaluateTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  ;

     virtual Real evaluateRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  ;

     virtual Real evaluateTimeDerivativeTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  ;

     virtual Real evaluateTimeDerivativeRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )   ;


    void setKernel(std::shared_ptr<kernel2B> kernel_) {_kernel=kernel_;}

    const auto & kernel() const {return *_kernel;}

    auto & kernel(){return *_kernel;}

    auto setA() const {return _setA;}
    auto setB() const  {return _setB;}
    
    void setSets(const std::array<int,2> & sets){
        _setA=sets[0];
        _setB=sets[0];
        
    }
    private:

    std::shared_ptr<kernel2B> _kernel;

    private:

    int _setA;
    int _setB;


};

struct twoBodySkippingEvaluationPolicy : public twoBodyEvaluationPolicy
{



     virtual Real evaluateTriangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    )  override ;
    
/*
    virtual Real evaluateRectangular(
        const configurations_t & configurations,
        const  std::array<int,2> & timeRange, 
        const std::array<int,2> & rangeA, 
        const std::array<int,2> & rangeB
    ) const override;
*/
    
};





class actionTwoBody : public action
{
    using range_t = std::array<int,2>;

    public:
    actionTwoBody() :
    _distanceMinimumConstraint(false),rMin(0)  {
        evaluationPolicy=std::make_shared<twoBodyEvaluationPolicy>();

    }

    void setKernel( std::shared_ptr<kernel2B> kernel_)
    {
        _kernel=kernel_;
        setGeometry(_kernel->geometry() );
        evaluationPolicy->setKernel(kernel_);
    }
    

    auto & kernel(){return _kernel;}
    const auto & kernel()const {return _kernel;}

    void setTwoBodyEvaluationPolicy(std::shared_ptr<twoBodyEvaluationPolicy> policy)
    {
        evaluationPolicy=policy;
    }

    void setSets(const std::array<int,2> & sets)
    {
        setA=sets[0];
        setB=sets[1];
        evaluationPolicy->setSets(sets);
    }


    Real evaluate(const configurations_t & configurations,const std::array<int,2> timeRange, int iParticle)
    {
        return evaluate(configurations,timeRange,{iParticle,iParticle});
    }


    Real evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange, const std::array<int,2> & particleRange);
    

    virtual Real evaluateTimeDerivative(const configurations_t & configurations,const std::array<int,2> & timeRange, const std::array<int,2> & particleRange) override;    

    virtual Real evaluate(const  configurations_t & configurations)
    {
        return evaluate(configurations,
                        range_t{0,configurations.nBeads()-1},
                        configurations.getGroups()[setA].range()
                        );
    }

    virtual Real evaluateTimeDerivative( const configurations_t & configurations) override
    {
        return evaluateTimeDerivative(configurations,
                        range_t{0,configurations.nBeads()-1},
                        configurations.getGroups()[setA].range()
                        );
    }
    

    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer);


    virtual void addGradient(const configurations_t &configurations,  Eigen::Tensor<Real,3> & gradientBuffer)
    {
        addGradient(configurations,
                        range_t{0,configurations.nBeads()-1},
                        configurations.getGroups()[setA].range(),
                        gradientBuffer
                        );
    }

    
    virtual bool checkConstraints(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange);
    



    void setMinimumDistance(Real rMin_) {
        _distanceMinimumConstraint=true;
        rMin=rMin_;
    };

    private:

    std::shared_ptr<kernel2B> _kernel;
    int setA;
    int setB;

    bool _distanceMinimumConstraint;

    std::shared_ptr<twoBodyEvaluationPolicy> evaluationPolicy;

    Real rMin;


};


}


#endif