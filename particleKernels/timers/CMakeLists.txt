include(../cmake/init.cmake)

add_executable(timers main.cpp)

if(OpenMP_CXX_FOUND)
target_link_libraries(timers PUBLIC OpenMP::OpenMP_CXX)
endif()
if (MPI_FOUND)
target_link_libraries(timers PUBLIC MPI::MPI_CXX)
endif()


target_include_directories(timers PUBLIC
                           "${PROJECT_SOURCE_DIR}/src"
			   )

target_link_libraries(timers PUBLIC particleKernels_lib )		
set_target_properties( timers PROPERTIES COMPILE_FLAGS ${CMAKE_CXX_FLAGS_ADDITIONAL} )
set_target_properties( timers PROPERTIES LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS_ADDITIONAL} )
