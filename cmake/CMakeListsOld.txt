cmake_minimum_required(VERSION 3.10)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if ( ${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_CXX_FLAGS_ADDITIONAL " -g -pg -Wfatal-errors")
  elseif (${CMAKE_BUILD_TYPE} MATCHES Release) 
  set(CMAKE_CXX_FLAGS_ADDITIONAL " -DNDEBUG -pg -Wfatal-errors")
  else()
  message(FATAL_ERROR "Unrecognized build type: " ${CMAKE_BUILD_TYPE}  )

endif()


# set the project name
project(qmc VERSION 3.0)

if(NOT Boost_INCLUDE_DIR)
  set(Boost_INCLUDE_DIR ${PROJECT_BINARY_DIR}/../external/boost )
endif()

find_package(MPI REQUIRED)
find_package(Boost REQUIRED)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)

include_directories(SYSTEM ${HDF5_INCLUDE_DIR})


set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})

include_directories(SYSTEM ${MPI_INCLUDE_PATH})
include_directories(SYSTEM ${Boost_INCLUDE_DIRS})


add_subdirectory(${PROJECT_BINARY_DIR}/../external/googletest)
add_subdirectory(${PROJECT_BINARY_DIR}/../external/eigen)

set(SOURCES  src/input.cpp src/geometry.cpp src/wavefunction/jastrows/jastrowSquareWell.cpp src/wavefunction/jastrows/jastrowSpline.cpp src/wavefunction/jastrows/jastrowHardSphereGauss.cpp src/wavefunction/jastrows/jastrowGaussGauss.cpp src/wavefunction/jastrows/jastrowPoschTeller.cpp src/wavefunction/jastrowWavefunctionThreeBody.cpp src/wavefunction/wavefunction.cpp src/wavefunction/productWavefunction.cpp src/tableDistances.cpp src/potential.cpp src/walkers.cpp src/initializer.cpp src/energy.cpp src/pairCorrelation.cpp src/trimerCorrelation.cpp src/estimators.cpp src/driver.cpp src/estimatorCollection.cpp src/tools.cpp src/moves/vmcMoves.cpp src/moves/dmcMoves.cpp src/dmcDriver.cpp src/vmcDriver.cpp src/metropolis.cpp src/branching.cpp src/orbitals.cpp src/slaters.cpp src/wavefunction/slaterDeterminant.cpp src/wavefunction/shell.cpp src/factory.cpp src/abstractFactory.cpp src/ptools.cpp src/timer.cpp src/estimatorsFactory.cpp src/correlationEstimator.cpp src/centerOfMassSquared.cpp src/hyperRadius.cpp src/wavefunction/jastrows/jastrow_delta.cpp src/wavefunction/jastrows/jastrowDipolar.cpp src/structureFactor.cpp src/parameters.cpp src/optimizationObservables.cpp src/geometryPMC.cpp
 pimc/action.cpp pimc/moves.cpp pimc/pimcConfigurations.cpp pimc/pimcObservables.cpp pimc/toolsPimc.cpp pimc/pimcDriver.cpp pimc/pimcPotentials.cpp pimc/hdf5IO.cpp
 
 )


add_library(main_lib STATIC ${SOURCES})
set_target_properties(main_lib PROPERTIES RULE_MESSAGES OFF)
set_target_properties( main_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( main_lib PROPERTIES LINK_FLAGS "-pg" )


target_include_directories(main_lib PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
			   )

target_include_directories(main_lib SYSTEM PUBLIC ${PROJECT_BINARY_DIR}/../external/eigen ${PROJECT_BINARY_DIR}/../external/json/single_include		           )
			 
target_link_libraries (main_lib eigen)
target_link_libraries (main_lib ${MPI_LIBRARIES} )
target_link_libraries( main_lib ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})


add_executable(main src/main.cpp)
target_link_libraries(main main_lib)
target_link_libraries(main )
set_target_properties( main PROPERTIES LINK_FLAGS "-pg" )
set_target_properties( main PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )



add_executable(pimc pimc/pimc.cpp)
target_include_directories(pimc PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           )
target_link_libraries(pimc PUBLIC  main_lib ${MPI_LIBRARIES})			 

add_executable(current_test src/current_test.cpp)
target_link_libraries(current_test main_lib ${MPI_LIBRARIES} )

set(TEST_SOURCES  tst/parameters-test.cpp tst/jastrow-test.cpp tst/wavefunction-test.cpp tst/dmc-test.cpp tst/ptools-test.cpp tst/pimc-test.cpp)
add_executable(tst tst/main.cpp ${TEST_SOURCES})
target_include_directories(tst PUBLIC
                           "${PROJECT_BINARY_DIR}/../src"
                           )
target_link_libraries(tst PUBLIC ${PROJECT_BINARY_DIR}/../external gtest eigen main_lib ${MPI_LIBRARIES})			 
set_target_properties( tst PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( tst PROPERTIES LINK_FLAGS "-pg" )


message("Compile Flags: " ${CMAKE_CXX_COMPILE_FLAGS} )