include(../cmake/init.cmake)

add_subdirectory("${PROJECT_SOURCE_DIR}/external/googletest" "external/googletest")
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

set(SOURCES  main.cpp verlet-test.cpp linkedCell.cpp twoBodyTest.cpp testUtils.cpp)


add_executable(tst ${SOURCES} )
target_include_directories(tst PUBLIC
                           "${PROJECT_SOURCE_DIR}/src"
			   )





if(OpenMP_CXX_FOUND)
target_link_libraries(tst PUBLIC OpenMP::OpenMP_CXX)
endif()
if (MPI_FOUND)
target_link_libraries(tst PUBLIC MPI::MPI_CXX)
endif()
target_link_libraries(tst PRIVATE gtest)
target_link_libraries(tst PUBLIC particleKernels_lib)


set_target_properties( particleKernels_lib PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( particleKernels_lib PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )