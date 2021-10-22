#include "testConfigurations.h"
#include "../pimc/toolsPimcTest.h"
#include "tools.h"
#include "../pimc/action.h"
#include "../pimc/pimcConfigurations.h"
#include "../pimc/moves.h"
#include "../pimc/pimcObservables.h"
#include "../pimc/hdf5IO.h"
#include "../pimc/nConnectedChains.h"
#include "../pimc/toolsPimcTest.h"

class twoBodyTest  : public configurationsTest
{
    public:

    void sample();

    virtual void accumulate();

    void SetUpCaoBernePropagator(Real radius);
    void SetRandomMinimumDistance(Real radius, const std::array<Real,DIMENSIONS> box ={TRUNCATE_D(1,1,1)});
    void SetUpFreeActionWithHardSphereConstraint(Real a);

    protected:

    int nTrials=100000;
    int nBlocks=100000;

    void setRandomMinimumDistance( Real radius);

};

class pbcTest  : public twoBodyTest
{
    public:

};

class harmonicTrapTest : public twoBodyTest
{
    public:
    void SetUpTwoBodyInteractionHarmonicInTrap_kernel()  ;
};
