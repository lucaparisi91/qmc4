#include "actionTwoBody.h"
#include "potentialKernelConstructor.h"

namespace pimc{


    class actionTwoBodyConstructor
    {
        public:

        using geometry_t = pimc::geometryPBC_PIMC;


        void setTimeStep( Real timeStep) { _timeStep = timeStep;
        primitiveConstructor.setTimeStep(timeStep); }



        template<class V_t>
        void registerPotential(const std::string & name)
        {
            primitiveConstructor.registerPotential<V_t>(name);

        }

        void setGeometry( const geometry_t & geo ) { _geo=geo;
        primitiveConstructor.setGeometry(geo); }


        auto create( const json_t & j )
        {
            auto kernel = primitiveConstructor.create(j);

            auto S = std::make_shared<actionTwoBody>() ;
            S->setKernel(kernel);

            auto setA = j["setA"].get<int>();
            auto setB = j["setB"].get<int>();

            S->setSets( { setA, setB});

            return S;

        }


        private:

        Real _timeStep;
        geometry_t _geo;

        primitiveApproximationTwoBodyKernelConstructor primitiveConstructor;

    };

}