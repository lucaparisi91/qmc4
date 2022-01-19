#include "exampleActions.h"



namespace pimc{

pimc::firstOrderAction createFreeCaoBerneAction(Real radius, Real timeStep , const geometry_t & geo,range_t groups, int nChains, int M)
{
    using propagator_t = pimc::caoBernePropagator;

    auto kernel =std::make_shared<propagator_t>(timeStep,radius);
    auto pairKernel = std::make_shared<pimc::pairProductKernel<propagator_t> >(kernel);
    auto buffer= std::make_shared<Eigen::Tensor<Real,3> >( nChains , getDimensions() , M + 1);

    pairKernel->setDistanceBuffer(buffer);

    auto S2B = std::make_shared<pimc::actionTwoBody>();
    S2B->setTimeStep(timeStep);
    S2B->setGeometry(geo);
    S2B->setKernel(pairKernel);
    pairKernel->setGeometry(geo);
    pairKernel->setTimeStep(timeStep);


    S2B->setSets({groups[0],groups[1]});
    S2B->setMinimumDistance(radius);
    S2B->setGeometry(geo);
    S2B->setTimeStep(timeStep);


    auto  sT= std::make_shared<pimc::kineticAction>(timeStep, nChains , M  , geo);


    std::shared_ptr<pimc::action> sEmpty= std::make_shared<pimc::nullPotentialAction>(timeStep  , geo);
    


    std::vector<std::shared_ptr<pimc::action> > Vs = {S2B};

    std::shared_ptr<pimc::action>  sV = std::make_shared<pimc::sumAction>(Vs);

    return pimc::firstOrderAction(sT,sV);


}

};