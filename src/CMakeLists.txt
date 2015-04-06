cmake_minimum_required (VERSION 2.6) 
project (HELLO) 

# set (CMAKE_Fortran_FLAGS_RELEASE "-O3")


SET ( SOURCES integrate.cpp lorenz63.c )



# set_directory_properties( PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES  "snap_shots;OUTPUT;parameters_output.txt;params.nml" )

# configure_file( deep_moist ${CMAKE_CURRENT_BINARY_DIR}/deep_moist COPYONLY)
add_library( lorenz63 SHARED ${SOURCES})
# add_executable ( lorenz63  ${SOURCES} )
target_link_libraries (lorenz63 m gsl gslcblas) 



