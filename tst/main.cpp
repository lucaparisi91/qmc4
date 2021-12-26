#include "gtest/gtest.h"
#include "ptools.h"
#include "test-tools.h"


int main(int argc, char **argv) {
  
  pTools::init(argc,argv);
  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  if (pTools::rank() != 0) {
    delete listeners.Release(listeners.default_result_printer());
}
  ::testing::InitGoogleTest(&argc, argv);

    std::string command_line_arg(argc == 2 ? argv[1] : "");

  ::testing::AddGlobalTestEnvironment(new testEnvironment(command_line_arg));
  

  int ret= RUN_ALL_TESTS();

  pTools::finalize();

  return ret;

}
