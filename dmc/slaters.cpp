#include "orbitals.h"
#include "slaters.h"


double getLog(double a,int & sign){if (a < 0 ) {sign*=-1;};return std::log(std::abs(a));}


std::complex<double> getLog(std::complex<double> a,int & sign){return std::log(a*(1.*sign));sign=1;}

void tableSlaters::add(int setA,orbitalSetBase * orbitalSet)
{
  
  if ( orbitalSet->isReal() )
    {
      
      int index = orbitalSetsReal.size();
      orbitalSetsReal.push_back(orbitalSet);
      indices1bReal[setA]=index;
      
      slaterMatricesReal.resize(slaterMatricesReal.size()+1);
      slaterMatricesRealInverse.resize(slaterMatricesRealInverse.size()+1);
      logDeterminantsReal.resize(logDeterminantsReal.size()+1);
      
      signsReal.resize(signsReal.size()+1);  
    }
  else
    {
      int index = orbitalSetsComplex.size();
      orbitalSetsComplex.push_back(orbitalSet);

      indices1bComplex[setA]=index;

      logDeterminantsComplex.resize(logDeterminantsComplex.size()+1);
      slaterMatricesComplex.resize( slaterMatricesComplex.size()+1);
      slaterMatricesComplexInverse.resize(slaterMatricesComplexInverse.size()+1);
      
      signsComplex.resize(signsComplex.size()+1);
  
    }
}

void tableSlaters::update(const tableSlaters::states_t & states)
{
  for ( const auto  element : indices1bReal ) // updates single set real slaters
	{
	  auto & matrix=slaterMatricesReal[element.second];
	  auto & matrixInverse=slaterMatricesRealInverse[element.second];
	  auto & logDeterminant= logDeterminantsReal[element.second];
	  auto & sign = signsReal[element.second];
	  
	  orbitalSetsReal[element.second]->storeEvaluate(states[element.first],matrix);
	  /* compute inverse */
	  ludReal.compute(matrix);
	  matrixInverse=ludReal.inverse();
	  
	  /* Compute determinant*/
	  logDeterminant=0;
	  auto & luMatrix=ludReal.matrixLU();
	  sign=ludReal.permutationP().determinant();
	  for(auto i=0;i<luMatrix.rows();i++)
	    {
	      logDeterminant+=getLog(luMatrix(i,i),sign );
	    }
	  
	  logDeterminantsReal[element.second]=logDeterminant;
	  

	}
  
    for ( const auto  element : indices1bComplex ) // updates single set real slaters
	{
	  auto & matrix=slaterMatricesComplex[element.second];
	  auto & matrixInverse=slaterMatricesComplexInverse[element.second];
	  auto & logDeterminant= logDeterminantsComplex[element.second];
	  auto & sign = signsComplex[element.second];
	  
	  orbitalSetsComplex[element.second]->storeEvaluate(states[element.first],matrix);
	  /* compute inverse */
	  ludComplex.compute(matrix);
	  matrixInverse=ludComplex.inverse();
	  
	  /* Compute determinant*/
	  logDeterminant=0;
	  auto & luMatrix=ludComplex.matrixLU();
	  sign=ludComplex.permutationP().determinant();
	  for(auto i=0;i<luMatrix.rows();i++)
	    {
	      logDeterminant+=getLog(luMatrix(i,i),sign );
	    }
	  
	  logDeterminantsComplex[element.second]=logDeterminant;
	  

	}

}
