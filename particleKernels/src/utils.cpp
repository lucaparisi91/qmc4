#include "utils.h"
#include <cmath>
#include <cassert>

namespace particleKernels
{

namespace utils
{

    Real restrictToBox(Real x, Real left, Real lBox, Real lBoxInverse)
    {

        x= x - std::round( ((x-left)*lBoxInverse - 0.5 ))*lBox;
        assert(x>=left);
        assert(x<=left+lBox);
        
        return x;
    };
    
    void restrictToBox(Real * positionsPBC, const Real * positionsOld, int iStart,int iEnd, int N, const std::array<Real,3> & left, const std::array<Real,3> & lBox )
    {
        std::array<Real,3> lBoxInverse{1./lBox[0],1./lBox[1],1./lBox[2]};
        for(int i=iStart;i<=iEnd;i++)
        {
            for(int d=0;d<3;d++)
            {
                positionsPBC[i+d*N]=restrictToBox(positionsOld[i+d*N],left[d],lBox[d],lBoxInverse[d]);
            }
        }
    }

    void initRandom(Real * particles,int dims,int iStart,int iEnd,size_t N,  std::default_random_engine & generator,
 Real xmin,Real xmax)
{

  std::uniform_real_distribution<double> distribution(xmin,xmax);


  for (int i=iStart;i<=iEnd;i++)
    {   
        for(int d=0;d<dims;d++)
        {
            particles[i+d*N]=distribution(generator);
        }
    }

}

};
}