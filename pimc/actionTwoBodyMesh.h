#ifndef ACTION_TWO_BODY_MESH_H
#define ACTION_TWO_BODY_MESH_H


#include "action.h"
#include "particleContainer/particleContainer.h"
#include "pairProductMeshKernel.h"
#include "pimcConfigurations.h"


namespace pimc
{
    struct actionTwoBodyMesh : public action
    {
        public:

        using geometry_t = geometryPBC_PIMC;

        actionTwoBodyMesh() {};

        virtual Real evaluate( const configurations_t & configurations, std::array<int,2> timeRange, int iParticle)
        {
            return _kernel->evaluate(configurations.dataTensor(),configurations.getAccelerationStructure().getLinkedCellList(),timeRange,iParticle);
        };

        virtual Real evaluate( const configurations_t & configurations)
        {
            const auto & group=configurations.getGroups()[0];

            return _kernel->evaluate(configurations.dataTensor(),configurations.getAccelerationStructure().getLinkedCellList(), {0,configurations.nBeads()-1},  {group.iStart,group.iEnd}   );
        };

        virtual Real evaluateTimeDerivative( const configurations_t & configurations)
        {
            const auto & group=configurations.getGroups()[0];

            return _kernel->evaluateTimeDerivative(configurations.dataTensor(),configurations.getAccelerationStructure().getLinkedCellList(), {0,configurations.nBeads()-1},  {group.iStart,group.iEnd}   );
        };

        void setKernel(std::shared_ptr<meshKernel> kernel_) {_kernel=kernel_; }

        auto & getKernel() { return *_kernel; }
        const auto &  getKernel() const { return *_kernel; }


        virtual void addGradient(const configurations_t & configurations,Eigen::Tensor<Real,3> & gradientBuffer) {
        
        const auto & group=configurations.getGroups()[0];
        
        _kernel->addForce(configurations.dataTensor(),configurations.getAccelerationStructure().getLinkedCellList(),{0,configurations.nBeads()-1},{group.iStart,group.iEnd},gradientBuffer);
        
        }

        

        private:

        std::shared_ptr<meshKernel> _kernel;

    };

}

#endif