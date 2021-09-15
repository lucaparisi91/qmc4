include(../cmake/init.cmake)

set(SOURCES  input.cpp geometry.cpp wavefunction/jastrows/jastrowSquareWell.cpp wavefunction/jastrows/jastrowSpline.cpp wavefunction/jastrows/jastrowHardSphereGauss.cpp wavefunction/jastrows/jastrowGaussGauss.cpp wavefunction/jastrows/jastrowPoschTeller.cpp wavefunction/jastrowWavefunctionThreeBody.cpp wavefunction/wavefunction.cpp wavefunction/productWavefunction.cpp tableDistances.cpp potential.cpp walkers.cpp initializer.cpp energy.cpp pairCorrelation.cpp trimerCorrelation.cpp estimators.cpp driver.cpp estimatorCollection.cpp tools.cpp moves/vmcMoves.cpp moves/dmcMoves.cpp dmcDriver.cpp vmcDriver.cpp metropolis.cpp branching.cpp orbitals.cpp slaters.cpp wavefunction/slaterDeterminant.cpp wavefunction/shell.cpp factory.cpp abstractFactory.cpp ptools.cpp timer.cpp estimatorsFactory.cpp correlationEstimator.cpp centerOfMassSquared.cpp hyperRadius.cpp wavefunction/jastrows/jastrow_delta.cpp wavefunction/jastrows/jastrowDipolar.cpp structureFactor.cpp parameters.cpp optimizationObservables.cpp
 )

 # -------  Set up the static dmc library
add_library(dmc_lib STATIC ${SOURCES})
include_qmc_external_directories(dmc_lib)
target_include_directories(dmc_lib PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
set_target_properties( dmc_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )


add_executable(dmc main.cpp)
target_include_directories(dmc PUBLIC
                           "${PROJECT_SOURCE_DIR}/dmc"
                           )
target_link_libraries(dmc PUBLIC dmc_lib)
set_target_properties( dmc PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} )
set_target_properties( dmc PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
link_qmc_external_libraries(dmc)
include_qmc_external_directories(dmc)
