#ifndef TRAITS_H
#define TRAITS_H
#include <vector>
#include <Eigen/Dense>
#include <random>
class geometry;
class productWavefunction;
#include "input.h"
#include "userDefinitions.h"



using real_t = double;
using geometry_t = geometry;
using state_t = Eigen::Array<real_t,Eigen::Dynamic,DIMENSIONS> ;
using states_t = std::vector<state_t> ;
class sumPotentials;

using difference_t = state_t;
using differences_t = std::vector<difference_t>;
using difference_complex_t = Eigen::Array<std::complex<real_t>,Eigen::Dynamic,DIMENSIONS> ;

using distance_t = Eigen::VectorXd;
using distances_t = std::vector<distance_t>;

using wavefunction_t = productWavefunction;

using randomGenerator_t = std::ranlux48;

//using json_t =   nlohmann::json;
using json_t = input;

using potential_t = sumPotentials;

#endif



namespace pimc
{
    class geometryPBC_PIMC;
    using Real = double;
    using geometry_t = geometryPBC_PIMC;
};