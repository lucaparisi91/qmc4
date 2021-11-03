#ifndef POTENTIAL_KERNEL_CONSTRUCTOR_H
#define POTENTIAL_KERNEL_CONSTRUCTOR_H



#include "potentialKernel.h"

namespace pimc
{

    namespace primitiveApproximationConstructors
    {

        struct constructorFromPotentialBase
        {
            virtual std::shared_ptr<kernel2B> create(const json_t & j) { throw std::runtime_error("Not implemented constructor from potential kernel"); };
        };

        template<class V_t> 
        struct constructorFromPotential : public          constructorFromPotentialBase
        {
            std::shared_ptr<kernel2B> create(const json_t & j)
            {
                auto V = std::make_shared<V_t>(j["potential"]);
                return std::make_shared<primitiveApproximationTwoBodyKernel<V_t> >(V);

            }
        };

    }


    class primitiveApproximationTwoBodyKernelConstructor
    {

        public:
        using geometry_t = geometryPBC_PIMC ;

        template<class V_t>
        void registerPotential( const std::string & name)
        {
            createorMap[name] =  primitiveApproximationConstructors::constructorFromPotential<V_t>();
        }

        void setTimeStep( Real timeStep){_timeStep=timeStep;}

        void setGeometry(const geometry_t & geo){_geo=geo;}

        auto create( const json_t & j )
        {
            auto key = j["potential"]["kind"];
            auto kernel= createorMap.at(key).create(j);
            kernel->setTimeStep(_timeStep);
            kernel->setGeometry(_geo);
            return kernel;
        }


        private:


        Real _timeStep;
        geometry_t _geo;


        std::map< std::string  ,  primitiveApproximationConstructors::constructorFromPotentialBase > createorMap;

    };

}

#endif