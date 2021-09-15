#include "gtest/gtest.h"
#include "testGlobals.h"


int main(int argc, char **argv) {
  
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
 /*  if (pTools::rank() != 0) {
    delete listeners.Release(listeners.default_result_printer());
} */
    ::testing::InitGoogleTest(&argc, argv);

    testGlobals::argc = argc;
    testGlobals::argv = argv;




  int ret= RUN_ALL_TESTS();


}
