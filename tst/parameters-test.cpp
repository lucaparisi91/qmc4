#include "gtest/gtest.h"
#include "wavefunction/jastrows/jastrow.h"
#include "wavefunction/productWavefunction.h"
#include <vector>
#include "traits.h"
#include "geometry.h"
#include "wavefunction/jastrowWavefunctionOneBody.h"
#include "wavefunction/jastrows/jastrow.h"
#include "parameters.h"
#include "wavefunction/jastrows/jastrowGaussian.h"


TEST(parametersTest,init)
{

  optimizationParameters parameters;
  bool status=parameters.addParameter("alpha",1);
  ASSERT_EQ(status,true);
  
  status=parameters.addParameter("alpha",1);
  ASSERT_EQ(status,false);
  
  optimizationParameter alpha("alpha",1,0);
  auto param = parameters.mapParameter(alpha,"alpha");
  
  
  gaussianJastrow J(1.0);
  
  
  std::vector<real_t> gradPs(1,0);
  
  J.addGradientParameter(1.,alpha,gradPs.begin(),gradPs.end() );
  
  ASSERT_EQ(gradPs[0],-1.);

  gradPs={0};
  
  J.addGradientParameter(1., param,gradPs);
  
  ASSERT_EQ(gradPs[0],-1.);
  
}

TEST(parametersTest,productWavefunction)
{
  int N=100;
	int D= getDimensions();
	
 	state_t particleData(N , D);
 	state_t gradient(N , D);
 	real_t alpha=1.;
 	
 	particleData.setRandom();
	
 	geometryPBC geo( 10., 10., 10.);
	

 	auto J=gaussianJastrow(alpha);
 	jastrowOneBodyWavefunction<gaussianJastrow> wave(J,geo,0);
	
 	productWavefunction waveT;
 	waveT.add(&wave);

}
