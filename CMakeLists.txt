cmake_minimum_required(VERSION 3.10)

# set the project name
project(qmc VERSION 3.0)

set(DIMENSIONS "1" CACHE STRING "Number of dimensions")
configure_file (
    "${PROJECT_SOURCE_DIR}/cmake/userDefinitions.h.in"
    "${PROJECT_SOURCE_DIR}/dmc/userDefinitions.h"
)

add_subdirectory(dmc)
enable_testing()
include(GoogleTest)
add_subdirectory(tst)
add_subdirectory(pimc)
add_subdirectory(particleKernels)