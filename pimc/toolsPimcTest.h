#ifndef TOOLS_PIMC_TEST
#define TOOLS_PIMC_TEST


#include "pimcConfigurations.h"

namespace pimc
{
    void testChain(pimc::pimcConfigurations & configurations, int iChain, int expectedHead, int expectedTail,const char * msg="\0");

    void testChainsJoined( pimcConfigurations & configurations, int iChainLeft, int iChainRight, const char * msg="\0");

#define TESTCHAIN(confs,iChain,head,tail) \
testChain(confs,iChain,head,tail, ("Chain test failed at file "  __FILE__  ":" + std::to_string(__LINE__ )).c_str()  )


#define TESTCHAINSJOINED(confs,iChain,iChainNext) \
testChainsJoined(confs,iChain,iChainNext, ("Chain test failed at file "  __FILE__  ":" + std::to_string(__LINE__ )).c_str()  )


};



#endif