#ifndef ACTION_TWO_BODY_MESH_H
#define ACTION_TWO_BODY_MESH_H


#include "action.h"
#include "particleContainer/particleContainer.h"
#include "pairProductMeshKernel.h"
#include "pimcConfigurations.h"

namespace pimc
{

    struct actionTwoBodyMesh
    {
        public:
        using geometry_t = geometryPBC_PIMC;

        actionTwoBodyMesh() {};

        
        virtual Real evaluate( const configurations_t & configurations, std::array<int,2> timeRange, int iParticle)
        {
            return _kernel->evaluate(configurations.dataTensor(),configurations.getAccelerationStructure().getLinkedCellList(),timeRange,iParticle);
        };


        void setKernel(std::shared_ptr<meshKernel> kernel_) {_kernel=kernel_; }

        auto & getKernel() { return *_kernel; }
        const auto &  getKernel() const { return *_kernel; }



        private:

        std::shared_ptr<meshKernel> _kernel;

    };

}

#endif