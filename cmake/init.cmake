

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(Boost REQUIRED)





function(include_qmc_external_directories target)



target_include_directories("${target}" SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories( ${target} PRIVATE ${PROJECT_SOURCE_DIR}/external/json/single_include)
target_include_directories(${target} PUBLIC ${PROJECT_SOURCE_DIR}/external)

endfunction()

function(link_qmc_external_libraries target)

string(REPLACE ":" " -L"  LIB_FLAGS "$ENV{LIBPATH}" )


#target_link_libraries (${target} PUBLIC eigen)
target_link_libraries(${target} PRIVATE hdf5)
target_link_libraries(${target} PUBLIC  particleKernels_lib )

if ( NOT (LIB_FLAGS STREQUAL "") )
  set_target_properties(${target} PROPERTIES LINK_FLAGS "-L ${LIB_FLAGS}" )
endif()



endfunction()


function(add_googletest)

add_subdirectory("${PROJECT_SOURCE_DIR}/external/googletest" "external/googletest")
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

endfunction()


add_compile_options(
  -Wfatal-errors
       $<$<CONFIG:RELEASE>:-O3>
       $<$<CONFIG:DEBUG>:-Og>
       $<$<CONFIG:DEBUG>:-g>
)

add_link_options(
       $<$<CONFIG:DEBUG>:-g>
)
add_compile_definitions(
        $<$<CONFIG:RELEASE>:NDEBUG>
        $<$<CONFIG:RELEASE>:BOOST_DISABLE_ASSERTS>
)