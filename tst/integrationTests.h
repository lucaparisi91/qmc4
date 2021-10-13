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

    protected:

    int nTrials=100000;
    int nBlocks=100000;



};

class pbcTest  : public twoBodyTest
{
    public:

    void sample();

    

};


