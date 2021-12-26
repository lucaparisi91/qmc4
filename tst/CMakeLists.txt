include(../cmake/init.cmake)

#add_googletest()

set(SOURCES_DMC   parameters-test.cpp jastrow-test.cpp wavefunction-test.cpp dmc-test.cpp ptools-test.cpp 

)


add_executable(tst_dmc main.cpp ${SOURCES_DMC})
target_include_directories(tst_dmc PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
target_link_libraries(tst_dmc PUBLIC dmc_lib)
target_link_libraries(tst_dmc PRIVATE gtest)
set_target_properties( tst_dmc PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
set_target_properties( tst_dmc PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
link_qmc_external_libraries(tst_dmc)
include_qmc_external_directories(tst_dmc)


set(SOURCES_PIMC PPA-test.cpp  pimc-test.cpp action-test.cpp testConfigurations.cpp integrationTests.cpp test-optimization.cpp acceleretorTest.cpp test-moves.cpp )

add_executable(tst_pimc main.cpp ${SOURCES_PIMC})
target_include_directories(tst_pimc PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
target_link_libraries(tst_pimc PUBLIC pimc_lib dmc_lib)
target_link_libraries(tst_pimc PRIVATE gtest)
set_target_properties( tst_pimc PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
set_target_properties( tst_pimc PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
link_qmc_external_libraries(tst_pimc)
include_qmc_external_directories(tst_pimc)
