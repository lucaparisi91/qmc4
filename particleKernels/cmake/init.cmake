
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)





find_package(OpenMP)


set(CMAKE_CXX_FLAGS "")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_compile_options(
  -Wfatal-errors
       $<$<CONFIG:RELEASE>:-O3>
       $<$<CONFIG:DEBUG>:-O0>
       $<$<CONFIG:DEBUG>:-g>
)


add_compile_definitions(
        $<$<CONFIG:RELEASE>:NDEBUG>
        $<$<CONFIG:RELEASE>:BOOST_DISABLE_ASSERTS>
)
