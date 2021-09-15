#include "jastrow.h"

class jastrowSpline : public  jastrow<jastrowSpline>
{
public:
  jastrowSpline(const json_t & j);
  
  jastrowSpline(Eigen::ArrayXd coefficients_,real_t stepSize_,real_t valueRight_,real_t derivativeRight=0);
  
  void evaluateDerivatives(real_t x, real_t & d0, real_t & d1, real_t & d2) ;

  real_t d0(real_t x);
  static std::string name() {return "bSpline"; }
  
private:
  
  void fillX(real_t x)
  {
    X[1]=x;
    X[2]=X[1]*x;
    X[3]=X[2]*x;
  }
  
  real_t stepSize;
  Eigen::ArrayXd coefficients;

  Eigen::Vector4d X;

  Eigen::Vector4d AX;
  Eigen::Vector4d AXd1;
  Eigen::Vector4d AXd2;
  
  const Eigen::Matrix4d A;
  const Eigen::Matrix4d Ad1;
  const Eigen::Matrix4d Ad2;
  
  real_t valueRight;
  real_t cutOff;
  real_t inverseStepSize;
  real_t derivativeRight;
  
};
