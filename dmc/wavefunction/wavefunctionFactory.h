#include "abstractFactory.h"

/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/

class wavefunction;

typedef wavefunction* (*wavefunctionCreatorFunc) ( const json_t & j ,const geometry_t & geo);

template<class wave_t>
wavefunction * createWavefunction(const json_t & j ,const geometry_t & geo)
{
  return new wave_t(j,geo);
}

class wavefunctionFactory : public abstractFactory<wavefunction,std::string, wavefunctionCreatorFunc>
{
public:
  using abstractFactory_t= abstractFactory<wavefunction,std::string, wavefunctionCreatorFunc>;

  
  template<class wave_t >
  void registerWavefunction()
  {
    registerType( wave_t::name()  , & (createWavefunction<wave_t> ) );
  }

  auto create(const json_t & j,const geometry_t & geo)
  {
    std::vector<wavefunction*> waves;

    for (const auto & waveJson : j )
      {
      
    

	std::string id = createId(waveJson);
	waves.push_back( abstractFactory_t::create(id,waveJson,geo));
      }
    
    return waves;
  }
  
 
  
};

