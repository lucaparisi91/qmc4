#ifndef KERNEL_CONSTRUCTOR_H
#define KERNEL_CONSTRUCTOR_H



#include "potentialKernel.h"
#include "pairProductKernel.h"

namespace pimc
{
    struct kernelConstructorBase
    {
        kernelConstructorBase() :
        _timeStep(0),_nBeads(0),_nMaxParticles(0) {}

        virtual std::shared_ptr<kernel2B> create(const json_t & j) { throw std::runtime_error("Not implemented constructor from potential kernel"); };

        void setTimeStep ( const Real & new_TimeStep_) { _timeStep=new_TimeStep_ ;};

        auto timeStep() { return _timeStep;}

        int nBeads() { return _nBeads;}
        int nMaxParticles() { return _nMaxParticles;}

        void setNBeads(int n) { _nBeads=n;}
        void setNMaxParticles(int n) { _nMaxParticles=n;}

        const auto & geometry(){return _geo;}
        void setGeometry(const geometry_t & geo){_geo=geo;}
        private:


        Real _timeStep;
        int _nBeads;
        int  _nMaxParticles;
        geometry_t _geo;
    };

    template<class V_t> 
    struct primitiveKernelConstructor : public          kernelConstructorBase
    {

         primitiveKernelConstructor() : kernelConstructorBase::kernelConstructorBase()
         {

         }

         void setPotential(std::shared_ptr<V_t> V) {_V=V;};

        auto create(  )
        {
            auto kernel= std::make_shared<primitiveApproximationTwoBodyKernel<V_t> >(_V);

            kernel->setGeometry( geometry());
            kernel->setTimeStep(  timeStep() );
            

            return kernel;

        }
        


        std::shared_ptr<kernel2B> create(const json_t & j)
        {
            auto V = std::make_shared<V_t>(j["potential"]);
            setPotential(V);

            return create();    
                      
        }

       

        private:

        std::shared_ptr<V_t> _V;
 
    };


    template<class G_t> 
    struct pairProductKernelConstructor : public   kernelConstructorBase {

        pairProductKernelConstructor() : kernelConstructorBase::kernelConstructorBase()
         {
             
         }

         
        void setGreenFunction( std::shared_ptr<G_t> G) { _G=G; }

        auto create()
        {
            auto kernel=std::make_shared< pairProductKernel<G_t> >(_G);       
             auto buffer = std::make_shared<Eigen::Tensor<Real,3> >(nMaxParticles()*nMaxParticles(),DIMENSIONS,nBeads()+1);
            buffer->setConstant(0.);
            kernel->setDistanceBuffer(buffer); 
            kernel->setGeometry( geometry());


            return kernel;

        }


        std::shared_ptr<kernel2B> create(const json_t & j)
        {
            json_t jG(j["greenFunction"]);
            jG["timeStep"]=timeStep();

            std::shared_ptr<G_t> G = std::make_shared<G_t>(jG);
            setGreenFunction(G);


            return create();


                  
        }

        private:

        std::shared_ptr<G_t> _G;


    };




    class kernelConstructor
    {

        public:
        using geometry_t = geometryPBC_PIMC ;

        kernelConstructor() :
        _timeStep(0) ,_nMaxParticles(0),_nBeads(0) {}

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
            
        }



        void setGeometry(const geometry_t & geo){_geo=geo;}


        auto create( const json_t & j )
        {
            checkInitialization();
            for ( auto it = createorMap.begin() ; it!= createorMap.end() ; it++ )
                {
                    it->second->setTimeStep(_timeStep);
                    it->second->setNBeads(_nBeads);
                    it->second->setNMaxParticles(_nMaxParticles);
                    it->second->setGeometry(_geo);

                }

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
            return kernel;
        }

        void setNBeads(int n){_nBeads=n;}
        void setNMaxParticles(int n){_nMaxParticles=n;}


        private:

        void checkInitialization(){
            if (_timeStep==0) throw std::runtime_error("timeStep non initialized in kernel constructor.");
            if (_nBeads==0) throw std::runtime_error("number of beads non initialized in kernel constructor.");
            if (_nMaxParticles==0) throw std::runtime_error("number of particles non initialized in kernel constructor.");

        }
        

        Real _timeStep;
        geometry_t _geo;
        int _nMaxParticles;
        int _nBeads;

        std::map< 
            std::string  ,  
            std::shared_ptr< kernelConstructorBase >    

        > createorMap;

    };

};


#endif