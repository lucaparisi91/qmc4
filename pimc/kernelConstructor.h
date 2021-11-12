#ifndef KERNEL_CONSTRUCTOR_H
#define KERNEL_CONSTRUCTOR_H



#include "potentialKernel.h"
#include "pairProductKernel.h"

namespace pimc
{
    struct kernelConstructorBase
    {
        virtual std::shared_ptr<kernel2B> create(const json_t & j) { throw std::runtime_error("Not implemented constructor from potential kernel"); };

        void setTimeStep ( const Real & new_TimeStep_) { _timeStep=new_TimeStep_ ;};

        auto timeStep() { return _timeStep;}
        private:

        Real _timeStep;
    };

    template<class V_t> 
    struct primitiveKernelConstructor : public          kernelConstructorBase
    {
        std::shared_ptr<kernel2B> create(const json_t & j)
        {
            auto V = std::make_shared<V_t>(j["potential"]);
            return std::make_shared<primitiveApproximationTwoBodyKernel<V_t> >(V);                
        }
    };


    template<class G_t> 
    struct pairProductKernelConstructor : public   kernelConstructorBase {
        std::shared_ptr<kernel2B> create(const json_t & j)
        {
            json_t jG(j["greenFunction"]);
            jG["timeStep"]=timeStep();

            auto G = std::make_shared<G_t>(jG);

            return std::make_shared< pairProductKernel<G_t> >(G);                
        }

    };




    class kernelConstructor
    {

        public:
        using geometry_t = geometryPBC_PIMC ;

        template<class V_t>
        void registerPotential( const std::string & name)
        {
            auto key = name + "__potential";

            createorMap[ key ] =  std::make_shared<
                primitiveKernelConstructor<V_t>
            > ();
        }

        template<class G_t>
        void registerGreenFunction( const std::string & name)
        {
            auto key = name + "__greenFunction";
            createorMap[ key ] =  std::make_shared<
                pairProductKernelConstructor<G_t>
            > ();
            createorMap[ key ]->setTimeStep(_timeStep);
        }

        void setTimeStep( Real timeStep){
            _timeStep=timeStep;
            for ( auto it = createorMap.begin() ; it!= createorMap.end() ; it++ )
                {
                    it->second->setTimeStep(timeStep);
                }
        }

        void setGeometry(const geometry_t & geo){_geo=geo;}


        auto create( const json_t & j )
        {
            std::string key="";
            if ( j.find("greenFunction") != j.end() )
            {
                key=j["greenFunction"]["kind"].get<std::string>() + "__greenFunction";

            }
            else if  ( j.find("potential") != j.end() )
            {
                key=j["potential"]["kind"].get<std::string>() + "__potential";

            }
            else
            {
                throw std::runtime_error("Green function of potential should be defined");
            }

            auto kernel= createorMap.at(key)->create(j);
            kernel->setTimeStep(_timeStep);
            kernel->setGeometry(_geo);
            return kernel;
        }


        private:


        Real _timeStep;
        geometry_t _geo;

        std::map< 
            std::string  ,  
            std::shared_ptr< kernelConstructorBase >    

        > createorMap;

    };

};


#endif