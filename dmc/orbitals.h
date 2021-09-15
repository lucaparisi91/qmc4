#include "traits.h"
#include <type_traits>
#include "qmcExceptions.h"
#include <complex>
#include <iomanip>

using namespace std::complex_literals;

class orbitalSetBase
{
public:
  using matrixReal_t = Eigen::MatrixXd;
  using matrixComplex_t = Eigen::Matrix<std::complex<real_t> , Eigen::Dynamic,Eigen::Dynamic>;
  
  virtual void storeEvaluate(const state_t & state,matrixReal_t & mat) const {throw missingImplementation("StoreEvaluate on a real matrix type");};
  virtual void storeEvaluate(const state_t & state,matrixComplex_t & mat) const { throw missingImplementation("StoreEvaluate on a complex matrix type");};
  
  virtual real_t energy() const =0;
  virtual bool isComplex() const = 0;
  virtual bool isReal() const = 0;
};


template<class orbital_t>
class orbitalSet : public orbitalSetBase
{
public:
  using value_t=typename orbital_t::value_t;
  using matrix_t = Eigen::Matrix<value_t,Eigen::Dynamic,Eigen::Dynamic>;
  
  
  virtual  bool isComplex() const {return  std::is_same<value_t, std::complex<real_t> >::value  ;};
  virtual bool isReal() const {return  std::is_same<value_t, real_t >::value  ;};
  
  orbitalSet() {};
  orbitalSet(const json_t & j);
  virtual void storeEvaluate(const state_t & state,matrix_t & matrix) const override; // store in a matrix all orbitals
  
  auto & getOrbitals() {return orbitals;}
  const auto & getOrbitals() const{return orbitals;}
  
  real_t energy() const 
  {
    real_t e=0;
    for ( auto & orbital : orbitals)
      {
	e+=orbital.energy();
      }
    return e;
  }

  static std::string name() {return orbital_t::name() ;}

private:
  std::vector<orbital_t> orbitals;
};

#if DIMENSIONS == 3


class sinOrbital
{
public:
  using value_t = real_t;
  static constexpr bool isComplex()  {return false;}
  
  sinOrbital(int nx,int ny,int nz,real_t lBox_);
  sinOrbital(int nx,int ny,int nz, const json_t & j);

  real_t operator()(real_t x,real_t y,real_t z) const { return std::sin(k[0]*x + k[1]*y + k[2]*z + delta) ; }
  
  inline void evaluateDerivatives(real_t x,real_t y, real_t z, real_t & dx, real_t & dy, real_t & dz ,real_t & laplacian ) const
  {
    real_t tmp = std::cos(k[0]*x + k[1]*y +k[2]*z + delta );
    real_t tmp1 = std::sin(k[0]*x + k[1]*y + k[2]*z + delta);
    dx=k[0]*tmp;
    dy=k[1]*tmp;
    dz=k[2]*tmp;
    laplacian=-(k[0]*k[0]+ k[1]*k[1]+ k[2]*k[2])*tmp1;
  }
  
  virtual real_t energy() const  {return 0.5*(k[0]*k[0] + k[1]*k[1] + k[2]*k[2] );}

  static std::string name() {return "sinCos" ;}
  
private:
  real_t delta;
  real_t lBox;
  std::vector<int> ns;
  std::vector<real_t> k;
  
};

class planeWave
{
public:
  using value_t = std::complex<real_t>;
  
  planeWave(int nx,int ny,int nz,real_t lBox_,real_t teta=0);
  planeWave(int nx,int ny, int nz, const json_t & json);
  static constexpr bool isComplex()  {return true;}
  
  std::complex<real_t> operator()(real_t x,real_t y,real_t z) const { return std::exp( std::complex<real_t>(0, ( k[0]*x + k[1]*y + k[2]*z)) ); }
  
  inline void evaluateDerivatives(real_t x,real_t y, real_t z, value_t & dx, value_t & dy, value_t & dz ,value_t & laplacian ) const
  {
    auto iImag = std::complex<real_t>(0,1);
    auto tmp = std::exp(std::complex<real_t>(0, (k[0]*x + k[1]*y +k[2]*z ) ) );
    dx=iImag*k[0]*tmp;
    dy=iImag*k[1]*tmp;
    dz=iImag*k[2]*tmp;
    laplacian= -(k[0]*k[0]+ k[1]*k[1]+ k[2]*k[2])*tmp;

  }
  
  virtual real_t energy() const  {return 0.5*(k[0]*k[0] + k[1]*k[1] + k[2]*k[2] );}
  
  static std::string name() {return "planeWave" ;}
  
private:
  real_t teta;
  real_t lBox;
  std::vector<int> ns;
  std::vector<real_t> k;
};







#include "wavefunction/shell.h"
template<class orbital_t, class ... Args >
void fillFermiSea(std::vector<orbital_t> & orbitals,int N,Args ... args)
{
  int n=0;
  shellStructure sh(N);
  
    for(int i=0;i<sh.nShells() and n<N;i++)
    {
      for(int j=0;j<sh[i].capacity() and n < N;j++)
	{
	  auto indices=sh[i][j];
	  orbitals.push_back(orbital_t(std::get<0>(indices),std::get<1>(indices),std::get<2>(indices), args ...  ));
	  n++;
	}
    }    
};



#include "wavefunction/shell.h"
template<class orbital_t, class ... Args >
void fillFermiSeaByEnergy(std::vector<orbital_t> & orbitals,int N,Args ... args)
{
  std::vector<orbital_t> orbitalsInShell;
  
  fillFermiSea(orbitalsInShell,2*N,args ...);
  std::sort(orbitalsInShell.begin(),orbitalsInShell.end(),[](orbital_t & a ,orbital_t & b){return a.energy() < b.energy();} );
  for (int i=0;i<N;i++)
   {
      orbitals.push_back(orbitalsInShell[i]);
    }
};


#endif
