#include "actionTwoBody.h"
#include "potentialKernelConstructor.h"
#include "actionConstructor.h"



namespace pimc{

    class actionTwoBodyConstructor : public actionConstructor
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


        std::shared_ptr<action> create( const json_t & j )
        {
            auto kernel = primitiveConstructor.create(j);

            auto S = std::make_shared<actionTwoBody>() ;
            S->setKernel(kernel);
            
            auto setA = j["groupA"].get<int>();
            auto setB = j["groupB"].get<int>();

            S->setSets( { setA, setB});
            S->setGeometry(_geo);
            S->setTimeStep(_timeStep);
            

            return S;

        }


        private:

        Real _timeStep;
        geometry_t _geo;

        primitiveApproximationTwoBodyKernelConstructor primitiveConstructor;

    };

}