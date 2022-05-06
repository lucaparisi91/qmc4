
#include "actionTwoBodyMesh.h"
#include "actionConstructor.h"
namespace pimc{


    class actionTwoBodyMeshConstructor : public actionConstructor
    {
        public:

        using geometry_t = pimc::geometryPBC_PIMC;

        actionTwoBodyMeshConstructor() : _timeStep(0),_nBeads(0),_nMaxParticles(0) {}
        
        void setTimeStep( Real timeStep) { _timeStep = timeStep;
         }

        void setGeometry( const geometry_t & geo ) { _geo=geo;}

        void setNMaxParticles( int n ) { _nMaxParticles=n; }

        void setNBeads(int n) { _nBeads=n;}


        std::shared_ptr<action> create( const json_t & j );

        private:

        Real _timeStep;
        geometry_t _geo;
        int _nMaxParticles;
        int _nBeads;


    };

}