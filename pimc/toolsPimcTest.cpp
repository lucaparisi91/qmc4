
#include "toolsPimcTest.h"
#include "gtest/gtest.h"

namespace pimc
{

void testChain(pimc::pimcConfigurations & configurations, int iChain, int expectedHead, int expectedTail,const char * msg)
{

    const auto chain = configurations.getChain(iChain);

    const auto & tags = configurations.getTags(); 

    ASSERT_EQ( chain.head , expectedHead   )<< std::endl <<  msg;
    ASSERT_EQ( chain.tail , expectedTail   )<< std::endl <<  msg;

    if (expectedHead > expectedTail)
    {
        for(int t=expectedTail+1;t<expectedHead;t++)
        {
            ASSERT_EQ(tags(iChain,t) , 1  ) << std::endl << "Failed at time " << t << std::endl <<  msg;
        }
        for(int t=0;t<=expectedTail;t++)
        {
            ASSERT_EQ(tags(iChain,t) , 0  )<< std::endl<< "Failed at time " << t << std::endl <<  msg;
        }
        for(int t=expectedHead;t<configurations.nBeads();t++)
        {
            ASSERT_EQ(tags(iChain,t) , 0  )<< std::endl<< "Failed at time " << t << std::endl <<  msg;
        }

    }
    else 
    {
        for(int t=expectedTail;t<=expectedHead;t++)
        {
            ASSERT_EQ(tags(iChain,t) , 0  )<< std::endl << "Failed at time " << t << std::endl<< msg;
        }
        for(int t=0;t<expectedTail;t++)
        {
            ASSERT_EQ(tags(iChain,t) , 1  )<< std::endl << "Failed at time " << t << std::endl<< msg;
        }
        for(int t=expectedHead+1;t<configurations.nBeads();t++)
        {
            ASSERT_EQ(tags(iChain,t) , 1  )<< std::endl << "Failed at time " << t << std::endl<< msg;
        }

    }

}


void testChainsJoined( pimcConfigurations & configurations, int iChainLeft, int iChainRight, const char * msg)
{
    const auto chainLeft = configurations.getChain(iChainLeft);
    const auto chainRight = configurations.getChain(iChainRight);

    // test if compatible with presence of heads or tails    
    ASSERT_FALSE(chainLeft.hasHead() ) << std::endl <<  msg;
    ASSERT_FALSE(chainRight.hasTail()  ) <<  std::endl <<  msg;

    // check that the link is bidirectional
    ASSERT_EQ(chainLeft.next , iChainRight)<<  std::endl <<  msg;
    ASSERT_EQ(chainRight.prev , iChainLeft)<<  std::endl <<  msg;


    // test that the two chains are part of the same group
    const auto & groupLeft = configurations.getGroupByChain(iChainLeft);
    const auto & groupRight = configurations.getGroupByChain(iChainRight);

    ASSERT_EQ(groupLeft.iStart , groupRight.iStart)<<  std::endl <<  msg;
    ASSERT_EQ(groupLeft.iEnd , groupRight.iEnd)<<  std::endl <<  msg;
}


}