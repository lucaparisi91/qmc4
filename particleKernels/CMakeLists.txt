cmake_minimum_required(VERSION 3.10)
# set the project name
project(particleKernels VERSION 1.0)

add_subdirectory(src)

enable_testing()
include(GoogleTest)
add_subdirectory(tst)

add_subdirectory(timers)
