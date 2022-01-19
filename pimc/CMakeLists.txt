include(../cmake/init.cmake)

set ( SOURCES geometryPMC.cpp
  action.cpp  moves.cpp  pimcConfigurations.cpp  pimcObservables.cpp pimcObservablesFactory.cpp toolsPimc.cpp  pimcDriver.cpp  pimcPotentials.cpp  hdf5IO.cpp toolsPimcTest.cpp
  forces.cpp propagators.cpp levyReconstructor.cpp  actionTwoBody.cpp nConnectedChains.cpp kernels.cpp accelerationStructure.cpp
  exampleActions.cpp) 

 # -------  Set up the static dmc library
add_library(pimc_lib STATIC ${SOURCES})
include_qmc_external_directories( pimc_lib)
target_include_directories( pimc_lib PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
set_target_properties(  pimc_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
target_link_libraries(pimc_lib PRIVATE gtest)


add_executable(pimc pimc.cpp)
target_include_directories(pimc PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
target_link_libraries(pimc PUBLIC stdc++fs  pimc_lib dmc_lib)



set_target_properties( pimc PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
target_link_libraries(pimc PRIVATE gtest )
set_target_properties( pimc PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
link_qmc_external_libraries(pimc)
include_qmc_external_directories(pimc)
