#include "action.h"
#include "actionConstructor.h"
namespace pimc
{

    namespace potentialActionOneBodyConstructors
    {
        struct potentialActionOneBodyConstructorBase
        {
            using geometry_t = geometryPBC_PIMC;
            void setTimeStep( Real timeStep) {_timeStep=timeStep; }
            void setGeometry(const geometry_t & geo) {_geo=geo; }
            virtual std::shared_ptr<action> create(const json_t & j){ throw std::runtime_error("Creation of potentialActionOneBody not implemented"); }
            protected:
            
            geometry_t _geo;
            Real _timeStep;

        };

        template<class V_t>
        struct potentialActionOneBodyConstructor : public potentialActionOneBodyConstructorBase
        {
            std::shared_ptr<action> create(const json_t & j)
            {
                auto set = j["group"].get<int>();

                V_t V(j["potential"]);


                return std::make_shared< potentialActionOneBody<V_t> >(_timeStep,V,_geo,set);

            }
        };

    }


    class potentialActionOneBodyConstructor : public actionConstructor
    {
        public:
        using geometry_t = geometryPBC_PIMC;

        void setTimeStep( Real timeStep) {_timeStep=timeStep; }
        void setGeometry(const geometry_t & geo) {_geo=geo; }


        virtual std::shared_ptr<action> create(const json_t & j) override
        {

            auto name = j["potential"]["kind"].get<std::string>();
            auto & constructor = keyToConstructor.at(name);

            constructor->setTimeStep(_timeStep);
            constructor->setGeometry(_geo);

            return constructor->create(j);

        }

        template<class V_t>
        void registerPotential(const std::string & name)
        {
            keyToConstructor[name]=std::make_shared<potentialActionOneBodyConstructors::potentialActionOneBodyConstructor<V_t> >();

        }

        private:

        Real _timeStep;
        geometry_t _geo;

        std::map<
            std::string,
            std::shared_ptr<potentialActionOneBodyConstructors::potentialActionOneBodyConstructorBase > 
                >
                keyToConstructor;

    };


};