#include "actionTwoBodyMeshConstructor.h"
#include "caoBerne.h"

namespace pimc 
{
     std::shared_ptr<action> actionTwoBodyMeshConstructor::create( const json_t & j )
        {

            auto S = std::make_shared<actionTwoBodyMesh>() ;

            S->setGeometry(_geo);
            S->setTimeStep(_timeStep);

            auto scatteringLength = j["greenFunction"]["a"].get<Real>();

            auto G = std::make_shared<pimc::caoBernePropagator>( _timeStep,scatteringLength);
            
            auto kernel=std::make_shared<pimc::pairProductMeshKernel<pimc::caoBernePropagator> >(G,_nMaxParticles, _nBeads);
            kernel->setGeometry(_geo);

            S->setKernel(kernel);

            return S;
        }
}