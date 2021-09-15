#include "pimcConfigurations.h"
#include <memory>


#include "potentialKernels.hpp"



namespace pimc
{
class potential
    {
        public:

        using return_t = Real;

        potential(int set) : _set(set) {}
        virtual return_t evaluate( const configurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange )=0;


        auto getSet() const {return _set;}

        private:

        int _set;
    };


template<class V_t>
class oneBodyIsotropicPotentialSecondOrder : public potential
{
    public:

    oneBodyIsotropicPotentialSecondOrder(V_t V, int set) : _V(V),potential(set) {}

    return_t evaluate( const configurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange 
    )
    {
        const auto & group= configurations.getGroups()[getSet() ];

        std::array<int,2> particleRangeClipped = {std::max(group.iStart,particleRange[0]),std::min(group.iEnd,particleRange[1])};

        if ( configurations.getEnsamble() == ensamble_t::grandCanonical )
        {

          
            return pimc::kernels::potentials::reduceOneBodyPotentialIsotropic(_V, configurations.dataTensor(), timeRange, particleRangeClipped, configurations.getTags() );
            
         

        }
        else if ( configurations.getEnsamble() == ensamble_t::canonical )
        {
                return pimc::kernels::potentials::reduceOneBodyPotentials2OrderIsotropic( _V,configurations.dataTensor(),timeRange, particleRangeClipped);
        }
        else
        {
            return 0.;
        }

        return 0;
    }

    private:

    V_t _V;
    
};

template<class V_t>
class oneBodyIsotropicPotentialChin : public potential
{
    public:

    oneBodyIsotropicPotentialChin(V_t V, int set,std::array<Real,2> vCoeff ) : _V(V),potential(set),_vCoeff(vCoeff) {}


    return_t evaluate( const configurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange 
    )
    {
        const auto & group= configurations.getGroups()[getSet() ];

        std::array<int,2> particleRangeClipped = {std::max(group.iStart,particleRange[0]),std::min(group.iEnd,particleRange[1])};

        if ( configurations.getEnsamble() == ensamble_t::grandCanonical )
        {

            return pimc::kernels::potentials::reduceOneBodyPotentialChin<V_t>(_V, configurations.dataTensor(), timeRange, particleRangeClipped, _vCoeff,configurations.getTags() );


        }
        else if ( configurations.getEnsamble() == ensamble_t::canonical )
        {
                return pimc::kernels::potentials::reduceOneBodyPotentialChin( _V,configurations.dataTensor(),timeRange, particleRangeClipped,_vCoeff);
        }
        else
        {
            return 0.;
        }

        return 0;
    }

    private:

    V_t _V;
    std::array<Real,2> _vCoeff;

};



class potentials
{
public: 
    potentials() {}

    void add(std::shared_ptr<potential> pot){_potentials.push_back(pot);}

    Real evaluate( const pimcConfigurations_t & configurations , const std::array<int,2> & timeRange, const  std::array<int,2> & particleRange )
    {
        Real sum=0;
        for (const auto & pot : _potentials)
        {
            sum+=pot->evaluate(configurations,timeRange,particleRange);
        }
        return sum;
    }


private:
    std::vector<std::shared_ptr<potential> > _potentials;
};

}