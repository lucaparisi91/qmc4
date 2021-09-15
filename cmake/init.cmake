

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if ( ${CMAKE_BUILD_TYPE} MATCHES Debug)
  set(CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL " -Wfatal-errors")
  set(CMAKE_CXX_LINK_FLAGS_ADDITIONAL "  ")

  
  elseif (${CMAKE_BUILD_TYPE} MATCHES Release) 
  set(CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL " -DNDEBUG -Wfatal-errors")
  set(CMAKE_CXX_LINK_FLAGS_ADDITIONAL " ")

  else()
  message(FATAL_ERROR "Unrecognized build type: " ${CMAKE_BUILD_TYPE}  )

endif()



find_package(MPI REQUIRED)
find_package(Boost REQUIRED)
find_package(HDF5 COMPONENTS C HL REQUIRED)


set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${CMAKE_CXX_COMPILE_FLAGS_ADDITIONAL} )
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )



function(include_qmc_external_directories target)

target_include_directories("${target}" SYSTEM PUBLIC ${Boost_INCLUDE_DIRS})
target_include_directories(${target} SYSTEM PUBLIC ${PROJECT_SOURCE_DIR}/external/eigen ${PROJECT_SOURCE_DIR}/external/json/single_include)
target_include_directories(${target} SYSTEM PUBLIC ${HDF5_INCLUDE_DIR})
target_include_directories(${target} PUBLIC ${PROJECT_SOURCE_DIR}/external)
target_include_directories(${target} SYSTEM PUBLIC ${MPI_CXX_INCLUDE_DIRS})

endfunction()

function(link_qmc_external_libraries target)

#target_link_libraries (${target} PUBLIC eigen)
target_link_libraries( ${target} PUBLIC ${HDF5_CXX_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(${target} PUBLIC MPI::MPI_CXX)
target_include_directories(${target} PUBLIC ${PROJECT_SOURCE_DIR}/external/eigen )
target_link_libraries(${target} PUBLIC  particleKernels_lib )


endfunction()


function(add_googletest)

add_subdirectory("${PROJECT_SOURCE_DIR}/external/googletest" "external/googletest")
mark_as_advanced(
    BUILD_GMOCK BUILD_GTEST BUILD_SHARED_LIBS
    gmock_build_tests gtest_build_samples gtest_build_tests
    gtest_disable_pthreads gtest_force_shared_crt gtest_hide_internal_symbols
)

endfunction()


message("Compile Flags: " ${CMAKE_CXX_COMPILE_FLAGS} )