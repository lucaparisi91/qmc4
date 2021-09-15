#ifndef FORCES_H
#define FORCES_H

#include "pimcConfigurations.h"


namespace pimc
{
    using pimcConfigurations_t = pimcConfigurations;


    class force
    {
        public:

        using tensor_t = pimcConfigurations_t::configurationsStorage_t ;

        virtual void add( const pimcConfigurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange, tensor_t & forces )=0;

    };

    template<class functor_t >
    class oneBodyIsotropicForce : public force
    {
        public:

        oneBodyIsotropicForce(functor_t V, int set) : _V(V),_set(set)
        {
            
        };


        virtual void add( const pimcConfigurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange, tensor_t & forces )
        {
            const auto & data = configurations.dataTensor();
            const auto & group = configurations.getGroups()[_set];


            std::array<int, 2 > intersectedParticleRange = { std::max(particleRange[0] , group.iStart) , std::min(particleRange[1],group.iEnd ) };


            for (int t=timeRange[0] ;t<= timeRange[1] ; t++  )
                for (int i=intersectedParticleRange[0];i<=intersectedParticleRange[1];i++ )
                    {
                        Real r2=0;
                        for (int d=0;d<getDimensions();d++)
                        {
                            r2+=data(i,d,t)*data(i,d,t);
                        }
                        Real r = std::sqrt(r2);

                        for (int d=0;d<getDimensions();d++)
                        {
                            forces(i,d,t)+=_V.radialDerivative(r)*data(i,d,t)/r;
                        }
                    }
        }

        private:

        functor_t _V;
        int _set;
    };


    class totalForce
    {
        public:

        totalForce(){}
        totalForce(int M , int nChains) { init(M,nChains); }
        void init(int M, int nChains);
        

        void evaluate(const pimcConfigurations_t & configurations,const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange);

        const auto & getForces() const {return forces; }


        const auto & getForceData() const {return *forceData; }
        auto & getForceData() {return *forceData; }


        using tensor_t = force::tensor_t;
        
        void add(std::shared_ptr<force> newForce){
            forces.push_back(newForce);
        }
        
        private:

        void setConstant(Real alpha,const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange);

        std::shared_ptr<tensor_t> forceData;
        std::vector<std::shared_ptr<force> > forces;
    };

}

#endif