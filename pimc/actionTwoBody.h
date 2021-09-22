#include "memory"
#include "kernels.h"
#include "action.h"

namespace pimc
{
//using configurations_t = pimcConfigurations;



class actionTwoBody : public action
{
    using range_t = std::array<int,2>;


    public:
    actionTwoBody()  {}

    void setKernel( std::shared_ptr<kernel2B> kernel_)
    {
        _kernel=kernel_;
        setGeometry(_kernel->geometry() );
    }


    void setSets(const std::array<int,2> & sets)
    {
        setA=sets[0];
        setB=sets[1];
    }


    Real evaluate(configurations_t & configurations,const std::array<int,2> timeRange, int iParticle)
    {
        return evaluate(configurations,timeRange,{iParticle,iParticle});
    }

    Real evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange, const std::array<int,2> & particleRange);

    Real evaluateTimeDerivative(const configurations_t & configurations,const std::array<int,2> & timeRange, const std::array<int,2> & particleRange);


    virtual Real evaluate( configurations_t & configurations)
    {
        return evaluate(configurations,
                        range_t{0,configurations.nBeads()-1},
                        configurations.getGroups()[setA].range()
                        );
    }

    virtual Real evaluateTimeDerivative( configurations_t & configurations)
    {
        return evaluateTimeDerivative(configurations,
                        range_t{0,configurations.nBeads()-1},
                        configurations.getGroups()[setA].range()
                        );
    }

    

    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer);


    private:

    std::shared_ptr<kernel2B> _kernel;
    int setA;
    int setB;



};

}