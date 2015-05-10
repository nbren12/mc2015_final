cmake_minimum_required (VERSION 2.6) 
project (HELLO) 

# Use CX11
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()


# Armadillo
# set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} /home/noah/usr/share/Armadillo/CMake/)
set(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH}
  /home/noah/usr/share/Armadillo/CMake
  /usr/local/share/Armadillo/CMake)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# set_directory_properties( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES  "snap_shots;OUTPUT;parameters_output.txt;params.nml" )
# configure_file( deep_moist ${CMAKE_CURRENT_BINARY_DIR}/deep_moist COPYONLY)

add_executable(param_search param_search.cpp equil.cpp integrate.cpp mcmc.cpp)
target_link_libraries (param_search m gsl gslcblas armadillo) 


