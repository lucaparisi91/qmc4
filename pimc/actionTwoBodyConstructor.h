#include "actionTwoBody.h"
#include "kernelConstructor.h"
#include "actionConstructor.h"


namespace pimc{

    class actionTwoBodyConstructor : public actionConstructor
    {
        public:


        using geometry_t = pimc::geometryPBC_PIMC;

        actionTwoBodyConstructor() : _timeStep(0),_nBeads(0),_nMaxParticles(0) {}

        void setTimeStep( Real timeStep) { _timeStep = timeStep;
         }


        template<class V_t>
        void registerPotential(const std::string & name)
        {
            kernelC.registerPotential<V_t>(name);
        }

        template<class G_t>
        void registerGreenFunction(const std::string & name)
        {
            kernelC.registerGreenFunction<G_t>(name);
        }

        void setGeometry( const geometry_t & geo ) { _geo=geo;}


        void setNMaxParticles( int n ) { _nMaxParticles=n; }

        void setNBeads(int n) { _nBeads=n;}
        

        std::shared_ptr<action> create( const json_t & j )
        {
            kernelC.setTimeStep(_timeStep);
            kernelC.setNBeads(_nBeads);
            kernelC.setGeometry(_geo);
            kernelC.setNMaxParticles(_nMaxParticles);


            auto kernel = kernelC.create(j);

            auto S = std::make_shared<actionTwoBody>() ;
            
            
            auto setA = j["groupA"].get<int>();
            auto setB = j["groupB"].get<int>();

            S->setSets( { setA, setB});
            S->setGeometry(_geo);
            S->setTimeStep(_timeStep);

            S->setKernel(kernel);
            kernel->setGeometry(_geo);
            kernel->setTimeStep(_timeStep);
            

            if ( j.find("minimumDistance") != j.end() )
            {
                Real radius = j["minimumDistance"].get<Real>();
                S->setMinimumDistance(radius);
            }
            

            return S;

        }


        private:

        Real _timeStep;
        geometry_t _geo;
        int _nMaxParticles;
        int _nBeads;

        kernelConstructor kernelC;

    };

}