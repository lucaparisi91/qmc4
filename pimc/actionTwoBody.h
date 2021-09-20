#include "memory"
#include "kernels.h"
#include "action.h"

namespace pimc
{
//using configurations_t = pimcConfigurations;


class actionTwoBody
{

    public:

    actionTwoBody() {}

    Real evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange, int iParticle);

    Real evaluate(const configurations_t & configurations,const std::array<int,2> & timeRange);
    
    virtual void addGradient(const configurations_t & pimcConfigurations,const std::array<int,2> & timeRange,const  std::array<int,2> & particleRange,  Eigen::Tensor<Real,3> & gradientBuffer){throw missingImplementation("Gradient not implemented for this action.");}



    
    private:

    


    std::shared_ptr<kernel2B> _kernel;
    int setA;
    int setB;



};

}