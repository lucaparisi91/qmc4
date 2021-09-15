#ifndef SLATERS_H
#define SLATERS_H


#include "traits.h"
#include <map>
#include <unordered_map>

class orbitalSetBase;

double getLog(double a,int & sign);
std::complex<double> getLog(std::complex<double> a,int & sign);

class tableSlaters
{
public:
  using matrixReal_t = Eigen::MatrixXd;
  using matrixComplex_t = Eigen::Matrix<std::complex<real_t>,Eigen::Dynamic,Eigen::Dynamic>;
  
  using states_t = ::states_t;
  
  tableSlaters(){};
  
  void update(const states_t & states); // update all registered distances
  
  void add(int setA ,orbitalSetBase * orbitalsSet); 
  
  //void add(int setA, int setB,orbitalSetBase * orbitalSet);
  
  const matrixReal_t & slaterMatrixReal(int setA) const { return slaterMatricesReal[indices1bReal.at(setA) ]; }
  const matrixReal_t & slaterMatrixRealInverse(int setA) const { return slaterMatricesRealInverse[indices1bReal.at(setA) ]; }

  const matrixComplex_t & slaterMatrixComplex(int setA) const { return slaterMatricesComplex[indices1bComplex.at(setA) ]; }
  const matrixComplex_t & slaterMatrixComplexInverse(int setA) const { return slaterMatricesComplexInverse[indices1bComplex.at(setA) ]; }


  
  //const slaterMatrix & slaterMatrix(int setA,int setB) const;
  
  real_t logDeterminantReal(int setA) const {return logDeterminantsReal[indices1bReal.at(setA)];}


  std::complex<real_t> logDeterminantComplex(int setA) const {return logDeterminantsComplex[indices1bComplex.at(setA)];}
  
  size_t size() {return orbitalSetsReal.size() + orbitalSetsComplex.size();}
  
private:
  
  const geometry_t *  geo;
  std::vector<orbitalSetBase* > orbitalSetsReal;
  std::vector<orbitalSetBase* > orbitalSetsComplex;
  
  
  std::vector<matrixReal_t> slaterMatricesReal;
  std::vector<matrixReal_t> slaterMatricesRealInverse;

  std::vector<matrixComplex_t> slaterMatricesComplex;
  std::vector<matrixComplex_t> slaterMatricesComplexInverse;
  
  
  std::unordered_map<int,int> indices1bReal;
  std::unordered_map<int,int> indices1bComplex;
  
  std::vector<real_t> logDeterminantsReal;
  std::vector<std::complex<real_t> > logDeterminantsComplex;
  
  std::vector<int> signsReal;
  std::vector<int> signsComplex;
  
  
  Eigen::PartialPivLU<matrixReal_t> ludReal;
  Eigen::PartialPivLU<matrixComplex_t> ludComplex;
  
  //std::map<std::pair<int,int> , int> indices2b;
};

#endif
