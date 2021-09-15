#include "abstractFactory.h"
#include "potentialTwoBody.h"
#include "potentialThreeBody.h"


/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/

class potential;

typedef potential* (*potentialCreatorFunc) ( const json_t & j ,const geometry_t & geo);

template<class pot_t>
potential * createPotential(const json_t & j ,const geometry_t & geo)
{
  return new pot_t(j,geo);
}


class potentialFactory : public abstractFactory<potential,std::string, potentialCreatorFunc>
{
public:
  using abstractFactory_t= abstractFactory<potential,std::string, potentialCreatorFunc>;

  
  template<class pot_t>
  void registerPotential()
  {
    registerType( pot_t::name()  , & (createPotential<pot_t> ) );
  }

  
  template<class functor_t>
  void registerFunctor()
  {
    registerPotential<potentialTwoBody<functor_t>  >();
    registerPotential<potentialThreeBodyUnDis<functor_t>  >();
    registerPotential<potentialThreeBodyDis<functor_t>  >();
  }
  
  auto create(const json_t & j,const geometry_t & geo)
  {
    std::vector<potential*> potentials;
    
    for (const auto & potJson : j )
      {
	std::string id = createId(potJson);
	potentials.push_back( abstractFactory_t::create(id,potJson,geo));
      }
    
    return potentials;
  }
  
 
  
};

