include(../cmake/init.cmake)

set(SOURCES  randomKernel.cpp cellNeighbourList.cpp timers.cpp
lattice.cpp particleData.cpp utils.cpp twoBodyPotential.cpp)

add_library(particleKernels_lib STATIC ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(particleKernels_lib PUBLIC OpenMP::OpenMP_CXX)
endif()
set_target_properties( particleKernels_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( particleKernels_lib PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
if (MPI_FOUND)
target_link_libraries(particleKernels_lib PUBLIC MPI::MPI_CXX)
endif()

add_executable(main main.cpp)

if(OpenMP_CXX_FOUND)
target_link_libraries(main PUBLIC OpenMP::OpenMP_CXX)
endif()
if (MPI_FOUND)
target_link_libraries(main PUBLIC MPI::MPI_CXX)
endif()

target_link_libraries(main PUBLIC particleKernels_lib)
set_target_properties( main PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( main PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )